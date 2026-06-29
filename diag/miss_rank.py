#!/usr/bin/env python3
"""
diag/miss_rank.py — where do the remaining misses actually sit?
==============================================================
For every current top-5 miss, deep-search the air-stripped apex spectrum to
top-30 and report the rank of the correct answer. This separates:
  rank 6-30  : answer IS in the pool -> recoverable (show top-10 / RI re-rank /
               soft tie-break) without new spectra
  not in 30  : answer absent -> needs better spectrum (deconv/apex) or library
and tags each miss type (alkane / oxygenated / long-chain / contaminant).
"""
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "eval"))

from config import (MZ_MIN, MZ_MAX, MZ_STEP, SMOOTH_WINDOW, MIN_SN,
                    MIN_PEAK_WIDTH_SCANS, ALKANE_RTS_DB_WAX)
from step1_parse import load_mzml_to_matrix
from step2_preprocess import estimate_column_bleed, subtract_column_bleed, preprocess_signal
from step3_peak_detect import detect_chromatographic_peaks
from step4_deconvolution import build_intensity_matrix
from step7_library_search import calc_ri_from_rt
from nist_engine import NISTSearchEngine, _is_alkane_like
from score import (names_match, load_ground_truth, load_candidates, score,
                   DEFAULT_GROUND_TRUTH, DEFAULT_CANDIDATES)
from bm_probe import find_peak, spec_apex_raw, spec_component, MZML


def tag(name):
    n = name.lower()
    if _is_alkane_like(name):
        return 'alkane'
    if 'hexadecen' in n or 'hexadecadien' in n:
        return 'longchain'
    if 'butyl' in n and 'phenol' in n:
        return 'contaminant'
    return 'oxygenated'


def main():
    user = load_ground_truth(DEFAULT_GROUND_TRUTH)
    peaks_json = load_candidates(str(DEFAULT_CANDIDATES))
    res = score(user, peaks_json)
    misses = [r['user'] for r in res if not r['top5'] and r['peak'] is not None]

    data = load_mzml_to_matrix(MZML)
    matrix, mz_bins = build_intensity_matrix(data['scan_list'], data['rt'],
                                             MZ_MIN, MZ_MAX, MZ_STEP)
    matrix = subtract_column_bleed(matrix, estimate_column_bleed(matrix), scale=0.5)
    pre = preprocess_signal(matrix.sum(axis=1), window_length=SMOOTH_WINDOW, max_half_window=100)
    peaks = detect_chromatographic_peaks(pre['corrected'], data['rt'], min_sn=MIN_SN,
                                         min_peak_width_scans=MIN_PEAK_WIDTH_SCANS,
                                         solvent_delay_min=4.0)
    engine = NISTSearchEngine()

    def rank(mz, inten, name, ri, depth=30):
        r = engine.search_raw_spectrum(mz, inten, top_n=depth, measured_ri=ri)
        for k, h in enumerate(r):
            if names_match(name, h['name']):
                return k + 1
        return None

    print("=" * 84)
    print(f"  {'RT':>6} {'compound':<34} {'type':>10} {'apex':>6} {'comp':>6}")
    print("=" * 84)
    buckets = {}
    for up in sorted(misses, key=lambda x: x['rt']):
        pk = find_peak(peaks, up['rt'])
        if pk is None:
            continue
        ri = calc_ri_from_rt(pk['apex_rt'], ALKANE_RTS_DB_WAX)
        ra = rank(*spec_apex_raw(data['scan_list'], pk), up['name'], ri)
        rc = rank(*spec_component(matrix, mz_bins, pk), up['name'], ri)
        t = tag(up['name'])
        best = min([x for x in (ra, rc) if x], default=None)
        key = 'in_top5' if (best and best <= 5) else \
              ('rank6_30' if best else 'absent')
        buckets.setdefault(key, []).append((up['name'], t, best))
        fa = f"#{ra}" if ra else "-"
        fc = f"#{rc}" if rc else "-"
        print(f"  {up['rt']:6.2f} {up['name'][:34]:<34} {t:>10} {fa:>6} {fc:>6}")
    print("=" * 84)
    for k in ('rank6_30', 'absent', 'in_top5'):
        items = buckets.get(k, [])
        print(f"  {k:>9}: {len(items)}")
        for nm, t, b in items:
            print(f"            {nm:<36} [{t}] best#{b}")
    print()
    print("  apex/comp = rank of correct answer in air-stripped top-30")
    print("  rank6_30 with component-spectrum win => deconvolution lever is real")


if __name__ == "__main__":
    main()
