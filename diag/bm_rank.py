#!/usr/bin/env python3
"""
diag/bm_rank.py — decisive Bm test: deep rank + characteristic-ion check.
=========================================================================
For each Bm compound, get NIST's TOP-50 candidates (apex + component spectra)
and report the deep rank of the correct answer:

  rank 1-5    : already a hit (shouldn't be here)
  rank 6-50   : RI re-ranking CAN rescue it  -> Phase 2 covers Bm too
  not in 50   : spectrum-data / ambiguity problem -> RI cannot help (M3)

Also prints the apex spectrum's top ions so missing characteristic ions
(mechanism M3) are visible.

Usage:  python diag/bm_rank.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "eval"))

from config import (MZ_MIN, MZ_MAX, MZ_STEP, SMOOTH_WINDOW, MIN_SN,
                    MIN_PEAK_WIDTH_SCANS, ALKANE_RTS_DB_WAX)
from step1_parse import load_mzml_to_matrix
from step2_preprocess import (preprocess_signal, estimate_column_bleed,
                              subtract_column_bleed)
from step3_peak_detect import detect_chromatographic_peaks
from step4_deconvolution import build_intensity_matrix
from step7_library_search import calc_ri_from_rt
from nist_engine import NISTSearchEngine
from score import names_match
from bm_probe import (BM_TARGETS, find_peak, spec_apex_raw, spec_component, MZML)


def deep_rank(engine, mz, inten, expected, ri, depth=50):
    res = engine.search_raw_spectrum(mz, inten, top_n=depth, measured_ri=ri)
    for k, r in enumerate(res):
        if names_match(expected, r['name']):
            return k + 1, res
    return None, res


def main():
    data = load_mzml_to_matrix(MZML)
    matrix, mz_bins = build_intensity_matrix(
        data['scan_list'], data['rt'], MZ_MIN, MZ_MAX, MZ_STEP)
    matrix = subtract_column_bleed(matrix, estimate_column_bleed(matrix), scale=0.5)
    pre = preprocess_signal(matrix.sum(axis=1), window_length=SMOOTH_WINDOW,
                            max_half_window=100)
    peaks = detect_chromatographic_peaks(
        pre['corrected'], data['rt'], min_sn=MIN_SN,
        min_peak_width_scans=MIN_PEAK_WIDTH_SCANS, solvent_delay_min=4.0)
    engine = NISTSearchEngine()

    print("=" * 78)
    print(f"  {'RT':>6} {'compound':<30} {'apex rank':>10} {'comp rank':>10}")
    print("=" * 78)
    verdict = {'rescue': [], 'absent': [], 'hit': []}
    for rt, name in BM_TARGETS:
        pk = find_peak(peaks, rt)
        if pk is None:
            continue
        ri = calc_ri_from_rt(pk['apex_rt'], ALKANE_RTS_DB_WAX)
        mzA, inA = spec_apex_raw(data['scan_list'], pk)
        mzC, inC = spec_component(matrix, mz_bins, pk)
        rA, _ = deep_rank(engine, mzA, inA, name, ri)
        rC, _ = deep_rank(engine, mzC, inC, name, ri)

        best = min([r for r in (rA, rC) if r], default=None)
        if best is None:
            verdict['absent'].append(name)
        elif best <= 5:
            verdict['hit'].append(name)
        else:
            verdict['rescue'].append((name, best))

        fa = f"#{rA}" if rA else "not in 50"
        fc = f"#{rC}" if rC else "not in 50"
        print(f"  {rt:6.2f} {name:<30} {fa:>10} {fc:>10}")

        # apex top ions (M3 visibility)
        order = np.argsort(inA)[::-1][:8]
        ions = "  ".join(f"{mzA[i]:.0f}:{inA[i]/inA.max()*100:.0f}" for i in order)
        print(f"         apex top ions  {ions}")

    print("=" * 78)
    print(f"  RESCUABLE by RI re-rank (rank 6-50): {len(verdict['rescue'])}")
    for nm, r in verdict['rescue']:
        print(f"       {nm}  (best rank #{r})")
    print(f"  ABSENT from top-50 (data/ambiguity, RI cannot help): "
          f"{len(verdict['absent'])}")
    for nm in verdict['absent']:
        print(f"       {nm}")
    print(f"  already hit: {len(verdict['hit'])}  {verdict['hit']}")


if __name__ == "__main__":
    main()
