#!/usr/bin/env python3
"""
diag/a_rank.py — Phase-2 gating test for A-class branched alkanes.
==================================================================
RI re-ranking can only rescue a compound if the CORRECT homolog is already
somewhere in NIST's candidate pool. For each A-class miss, fetch NIST's
top-50 (air-stripped apex spectrum) and report the deep rank of the exact
correct compound. Also report whether that compound has a DB-WAX literature
RI (if not, a PREDICTED RI is required -> the CNN model).

  rank 6-50 : RI re-rank can promote it  -> Phase 2 viable
  not in 50 : EI cannot place it in the pool -> RI cannot help

Usage:  python diag/a_rank.py
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
from score import names_match
from bm_probe import find_peak, spec_apex_raw, MZML

A_TARGETS = [
    (19.19, "2,4,6-Trimethyldecane"),
    (19.83, "2,6,10-Trimethyldodecane"),
    (21.01, "2,6,10,14-Tetramethylhexadecane"),
    (22.98, "2,6,10-Trimethylpentadecane"),
    (25.45, "2,6,10-Trimethyltetradecane"),
    (26.28, "2,6,10,14-Tetramethylhexadecane"),
    (26.66, "2,6,10-Trimethyldodecane"),
    (28.38, "2,6,11-Trimethyldodecane"),
    (30.54, "2,6,10,14-Tetramethylpentadecane"),
    (32.04, "2,6,10-Trimethylhexadecane"),
    (39.76, "2,6,10-Trimethyltridecane"),
]


def main():
    data = load_mzml_to_matrix(MZML)
    matrix, mz_bins = build_intensity_matrix(
        data['scan_list'], data['rt'], MZ_MIN, MZ_MAX, MZ_STEP)
    matrix = subtract_column_bleed(matrix, estimate_column_bleed(matrix), scale=0.5)
    pre = preprocess_signal(matrix.sum(axis=1), window_length=SMOOTH_WINDOW, max_half_window=100)
    peaks = detect_chromatographic_peaks(pre['corrected'], data['rt'], min_sn=MIN_SN,
                                         min_peak_width_scans=MIN_PEAK_WIDTH_SCANS,
                                         solvent_delay_min=4.0)
    engine = NISTSearchEngine()
    # has DB-WAX literature RI for the exact compound?
    wax = engine.ri_wax

    print("=" * 82)
    print(f"  {'RT':>6} {'compound':<32} {'rank':>9} {'litRI?':>7}  pool examples")
    print("=" * 82)
    rescuable = 0
    for rt, name in A_TARGETS:
        pk = find_peak(peaks, rt)
        if pk is None:
            print(f"  {rt:6.2f} {name:<32} [no peak]")
            continue
        ri = calc_ri_from_rt(pk['apex_rt'], ALKANE_RTS_DB_WAX)
        mz, inten = spec_apex_raw(data['scan_list'], pk)
        res = engine.search_raw_spectrum(mz, inten, top_n=50, measured_ri=ri)
        rank = None
        for k, r in enumerate(res):
            if names_match(name, r['name']):
                rank = k + 1
                break
        litri = 'Y' if name.lower() in wax else '-'
        if rank and 6 <= rank <= 50:
            rescuable += 1
        # show a few branched-alkane candidates already in the pool
        pool = [r['name'] for r in res[:50]
                if 'methyl' in r['name'].lower() and 'ane' in r['name'].lower()][:3]
        rs = f"#{rank}" if rank else "not in 50"
        print(f"  {rt:6.2f} {name:<32} {rs:>9} {litri:>7}  {'; '.join(pool)[:30]}")
    print("=" * 82)
    print(f"  Rescuable by RI re-rank (exact compound at rank 6-50): {rescuable}/{len(A_TARGETS)}")
    print(f"  litRI? = exact compound present in DB-WAX literature table "
          f"(else CNN-predicted RI needed)")


if __name__ == "__main__":
    main()
