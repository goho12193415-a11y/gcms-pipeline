#!/usr/bin/env python3
"""
diag/bgsub_test.py — does Xcalibur-style background subtraction recover misses?
==============================================================================
The native RAW spectrum is just as air/background-dominated as the mzML one, so
reading the RAW natively does not help by itself. What the analyst does in
Xcalibur/Chromeleon is BACKGROUND SUBTRACTION: pick a nearby baseline region,
average it, subtract it from the peak's apex spectrum. This removes air, column
bleed, and co-eluting baseline at once. Test whether that recovers the misses
beyond the current air-ion stripping.
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
from nist_engine import NISTSearchEngine
from score import load_ground_truth, load_candidates, score, DEFAULT_GROUND_TRUTH, DEFAULT_CANDIDATES
from bm_probe import find_peak, MZML


def bg_subtract(matrix, peak, width=10, gap=2):
    """Apex spectrum minus a local baseline spectrum taken just before the peak
    rises (a valley region), like Xcalibur's background subtraction."""
    a, s = peak['apex_idx'], peak['start_idx']
    lo = max(0, s - gap - width)
    hi = max(lo + 1, s - gap)
    bg = matrix[lo:hi, :].mean(axis=0) if hi > lo else np.zeros(matrix.shape[1])
    spec = matrix[a, :] - bg
    spec[spec < 0] = 0
    return spec


def main():
    user = load_ground_truth(DEFAULT_GROUND_TRUTH)
    res0 = score(user, load_candidates(str(DEFAULT_CANDIDATES)), top_n=5)
    misses = [r['user'] for r in res0 if not r['top5'] and r['peak'] is not None]

    data = load_mzml_to_matrix(MZML)
    matrix, mz_bins = build_intensity_matrix(data['scan_list'], data['rt'],
                                             MZ_MIN, MZ_MAX, MZ_STEP)
    matrix = subtract_column_bleed(matrix, estimate_column_bleed(matrix), scale=0.5)
    pre = preprocess_signal(matrix.sum(axis=1), window_length=SMOOTH_WINDOW, max_half_window=100)
    peaks = detect_chromatographic_peaks(pre['corrected'], data['rt'], min_sn=MIN_SN,
                                         min_peak_width_scans=MIN_PEAK_WIDTH_SCANS,
                                         solvent_delay_min=4.0)
    eng = NISTSearchEngine()
    from score import names_match

    def topk(spec, name, ri, k=5):
        mask = spec > 0
        if not mask.any():
            return None
        r = eng.search_raw_spectrum(mz_bins[mask], spec[mask], top_n=k, measured_ri=ri)
        for i, h in enumerate(r):
            if names_match(name, h['name']):
                return i + 1
        return (0, r[0]['name'] if r else '-')

    print("=" * 74)
    print(f"  {'RT':>6} {'compound':<34} {'bg-sub top-5':>12}")
    print("=" * 74)
    gained = 0
    for up in sorted(misses, key=lambda x: x['rt']):
        pk = find_peak(peaks, up['rt'])
        if pk is None:
            continue
        ri = calc_ri_from_rt(pk['apex_rt'], ALKANE_RTS_DB_WAX)
        spec = bg_subtract(matrix, pk)
        r = topk(spec, up['name'], ri)
        hit = isinstance(r, int)
        if hit:
            gained += 1
        tag = f"#{r}" if hit else "-"
        print(f"  {up['rt']:6.2f} {up['name'][:34]:<34} {tag:>12}")
    print("=" * 74)
    print(f"  recovered by background subtraction: {gained}/{len(misses)}")


if __name__ == "__main__":
    main()
