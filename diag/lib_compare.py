#!/usr/bin/env python3
"""
diag/lib_compare.py — score the 49 against an ARBITRARY NIST library path.
=========================================================================
Used to test whether the NIST26-EI-DEMO library (the one actually used to make
the manual IDs) matches better than the pipeline's NIST 2014 mainlib.

    python diag/lib_compare.py "C:\\NIST26-EI-DEMO\\MSSEARCH\\mainlib"
"""
import sys, json
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
from score import (load_ground_truth, score, DEFAULT_GROUND_TRUTH)
from bm_probe import MZML


def main():
    lib_path = sys.argv[1] if len(sys.argv) > 1 else None
    label = lib_path or "default (NIST2014 mainlib)"
    print(f"Library under test: {label}")

    data = load_mzml_to_matrix(MZML)
    matrix, mz_bins = build_intensity_matrix(data['scan_list'], data['rt'],
                                             MZ_MIN, MZ_MAX, MZ_STEP)
    matrix = subtract_column_bleed(matrix, estimate_column_bleed(matrix), scale=0.5)
    pre = preprocess_signal(matrix.sum(axis=1), window_length=SMOOTH_WINDOW, max_half_window=100)
    peaks = detect_chromatographic_peaks(pre['corrected'], data['rt'], min_sn=MIN_SN,
                                         min_peak_width_scans=MIN_PEAK_WIDTH_SCANS,
                                         solvent_delay_min=4.0)
    for p in peaks:
        p['ri_measured'] = calc_ri_from_rt(p['apex_rt'], ALKANE_RTS_DB_WAX)

    engine = NISTSearchEngine(lib_path=lib_path)
    matches = engine.search_peaks_raw(peaks, data['scan_list'], top_n=5, verbose=False)

    cand = [{'rt': round(p['apex_rt'], 3),
             'candidates': [{'name': m.get('name', ''), 'rmf': m.get('rmf', 0)}
                            for m in (matches[i] or [])[:5]]}
            for i, p in enumerate(peaks)]

    user = load_ground_truth(DEFAULT_GROUND_TRUTH)
    res = score(user, cand)
    n = len(res)
    t1 = sum(1 for r in res if r['top1'])
    t5 = sum(1 for r in res if r['top5'])
    print("=" * 56)
    print(f"  top-1: {t1}/{n} = {t1*100/n:.0f}%")
    print(f"  top-5: {t5}/{n} = {t5*100/n:.0f}%")
    print("=" * 56)
    print("  misses:")
    for r in sorted((x for x in res if not x['top5']), key=lambda x: x['user']['rt']):
        top = r['peak']['candidates'][0]['name'] if (r['peak'] and r['peak']['candidates']) else '-'
        print(f"    RT={r['user']['rt']:5.2f} {r['user']['name'][:32]:<32} -> {top[:34]}")


if __name__ == "__main__":
    main()
