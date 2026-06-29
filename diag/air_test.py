#!/usr/bin/env python3
"""Decisive test: does stripping air ions rescue the Bm class?"""
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
from bm_probe import BM_TARGETS, find_peak, spec_apex_raw, MZML

AIR = {17, 18, 28, 32, 40, 44}          # H2O, N2/CO, O2, Ar, CO2
LOW_CUTOFF = 41                          # alt: drop everything below m/z 41


def strip(mz, inten, mode):
    mz = np.asarray(mz, float); inten = np.asarray(inten, float)
    if mode == 'air':
        keep = np.array([round(m) not in AIR for m in mz])
    elif mode == 'low':
        keep = mz >= LOW_CUTOFF
    else:
        keep = np.ones(len(mz), bool)
    return mz[keep], inten[keep]


def rank_of(engine, mz, inten, expected, ri, depth=10):
    res = engine.search_raw_spectrum(mz, inten, top_n=depth, measured_ri=ri)
    for k, r in enumerate(res):
        if names_match(expected, r['name']):
            return k + 1, (res[0]['name'] if res else '-'), (res[0]['rmf'] if res else 0)
    return None, (res[0]['name'] if res else '-'), (res[0]['rmf'] if res else 0)


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
    print("=" * 76)
    print(f"  {'RT':>6} {'compound':<30} {'raw':>6} {'-air':>6} {'-<41':>6}")
    print("=" * 76)
    score = {'raw': 0, 'air': 0, 'low': 0}
    for rt, name in BM_TARGETS:
        pk = find_peak(peaks, rt)
        if pk is None:
            continue
        ri = calc_ri_from_rt(pk['apex_rt'], ALKANE_RTS_DB_WAX)
        mz, inten = spec_apex_raw(data['scan_list'], pk)
        r0, t0, m0 = rank_of(engine, *strip(mz, inten, 'none'), name, ri)
        ra, ta, ma = rank_of(engine, *strip(mz, inten, 'air'), name, ri)
        rl, tl, ml = rank_of(engine, *strip(mz, inten, 'low'), name, ri)
        for key, r in (('raw', r0), ('air', ra), ('low', rl)):
            if r and r <= 5:
                score[key] += 1
        f = lambda r: (f"#{r}" if r else "-")
        print(f"  {rt:6.2f} {name:<30} {f(r0):>6} {f(ra):>6} {f(rl):>6}")
        print(f"         -air top-1: {ta[:48]:<48} RMF={ma}")
    print("=" * 76)
    n = len(BM_TARGETS)
    print(f"  in top-5:  raw={score['raw']}/{n}   strip-air={score['air']}/{n}   "
          f"drop<41={score['low']}/{n}")


if __name__ == "__main__":
    main()
