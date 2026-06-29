#!/usr/bin/env python3
"""Test whether adding NIST replib recovers any of the 26 misses."""
import sys, os, tempfile
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "eval"))

from config import (MZ_MIN, MZ_MAX, MZ_STEP, SMOOTH_WINDOW, MIN_SN,
                    MIN_PEAK_WIDTH_SCANS, ALKANE_RTS_DB_WAX, AIR_IONS)
from step1_parse import load_mzml_to_matrix
from step2_preprocess import estimate_column_bleed, subtract_column_bleed, preprocess_signal
from step3_peak_detect import detect_chromatographic_peaks
from step4_deconvolution import build_intensity_matrix
from step7_library_search import calc_ri_from_rt
from score import names_match, load_ground_truth, load_candidates, score, DEFAULT_CANDIDATES
from bm_probe import find_peak, spec_apex_raw, MZML

import pyms_nist_search
from pyms.Spectrum import MassSpectrum

REPLIB = r"C:\Users\go ho\Desktop\MSSEARCH\replib"
AIR = frozenset(AIR_IONS)


def search(engine, mz, inten, name, depth=10):
    mz = np.asarray(mz, float); inten = np.asarray(inten, float)
    keep = np.array([round(m) not in AIR for m in mz])
    mz, inten = mz[keep], inten[keep]
    if len(mz) == 0 or inten.max() <= 0:
        return None
    order = np.argsort(mz)
    mz, inten = mz[order], (inten[order] / inten.max() * 999).astype(int)
    try:
        hits = engine.full_search_with_ref_data(MassSpectrum(mz.tolist(), inten.tolist()))
    except Exception:
        return None
    for k, h in enumerate(hits[:depth]):
        if isinstance(h, tuple) and len(h) >= 2:
            if names_match(name, getattr(h[1], 'name', '') or ''):
                return k + 1
    return None


def main():
    user = load_ground_truth(r"C:\Users\go ho\Desktop\桌面\xt_filled.xlsx")
    res = score(user, load_candidates(str(DEFAULT_CANDIDATES)))
    misses = [r['user'] for r in res if not r['top5'] and r['peak'] is not None]

    data = load_mzml_to_matrix(MZML)
    matrix, mz_bins = build_intensity_matrix(data['scan_list'], data['rt'],
                                             MZ_MIN, MZ_MAX, MZ_STEP)
    matrix = subtract_column_bleed(matrix, estimate_column_bleed(matrix), scale=0.5)
    pre = preprocess_signal(matrix.sum(axis=1), window_length=SMOOTH_WINDOW, max_half_window=100)
    peaks = detect_chromatographic_peaks(pre['corrected'], data['rt'], min_sn=MIN_SN,
                                         min_peak_width_scans=MIN_PEAK_WIDTH_SCANS,
                                         solvent_delay_min=4.0)

    wd = os.path.join(tempfile.gettempdir(), "replib_work"); os.makedirs(wd, exist_ok=True)
    print("Loading replib engine...")
    eng = pyms_nist_search.Engine(REPLIB, pyms_nist_search.NISTMS_REP_LIB, wd)

    print("=" * 70)
    print(f"  {'RT':>6} {'compound':<36} {'replib top-10':>13}")
    print("=" * 70)
    recovered5 = 0
    for up in sorted(misses, key=lambda x: x['rt']):
        pk = find_peak(peaks, up['rt'])
        if pk is None:
            continue
        mz, inten = spec_apex_raw(data['scan_list'], pk)
        r = search(eng, mz, inten, up['name'])
        if r and r <= 5:
            recovered5 += 1
        print(f"  {up['rt']:6.2f} {up['name'][:36]:<36} {('#'+str(r) if r else '-'):>13}")
    print("=" * 70)
    print(f"  replib would add to top-5: {recovered5}/{len(misses)} misses")


if __name__ == "__main__":
    main()
