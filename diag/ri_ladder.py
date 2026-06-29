#!/usr/bin/env python3
"""
diag/ri_ladder.py — extract the n-alkane RI ladder from the zgwt standard.
==========================================================================
The zgwt run is the C10-C27 n-alkane standard. Its major peaks ARE the
n-alkanes, eluting in carbon order. We detect them, confirm each is a
straight-chain alkane (base peak m/z 57, series 43/57/71/85...), read the
accurate apex RT, and assign carbon number via the molecular ion M+ = 14n+2
where visible (<= m/z 300 cutoff) and RT order otherwise.

Output: an accurate [(RT, RI=100*n), ...] ladder for calc_ri_from_rt.
"""
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import (MZ_MIN, MZ_MAX, MZ_STEP, SMOOTH_WINDOW, MIN_SN, MIN_PEAK_WIDTH_SCANS)
from step1_parse import load_sample
from step2_preprocess import estimate_column_bleed, subtract_column_bleed, preprocess_signal
from step3_peak_detect import detect_chromatographic_peaks
from step4_deconvolution import build_intensity_matrix

ZGWT = r"C:\Users\go ho\Desktop\验证\zgwt.raw"


def is_alkane(mz, inten):
    """Straight/branched alkane: base peak in {43,57,71}, strong CnH2n+1 series."""
    if len(mz) == 0:
        return False
    base = int(round(mz[np.argmax(inten)]))
    series = [43, 57, 71, 85, 99]
    mzr = np.round(mz).astype(int)
    hits = sum(1 for s in series if s in mzr)
    return base in (43, 57, 71, 85) and hits >= 4


def main():
    data = load_sample(ZGWT)
    matrix, mz_bins = build_intensity_matrix(data['scan_list'], data['rt'],
                                             MZ_MIN, MZ_MAX, MZ_STEP)
    matrix = subtract_column_bleed(matrix, estimate_column_bleed(matrix), scale=0.5)
    pre = preprocess_signal(matrix.sum(axis=1), window_length=SMOOTH_WINDOW, max_half_window=100)
    peaks = detect_chromatographic_peaks(pre['corrected'], data['rt'], min_sn=MIN_SN,
                                         min_peak_width_scans=MIN_PEAK_WIDTH_SCANS,
                                         solvent_delay_min=4.0)
    print(f"{len(peaks)} peaks in zgwt")

    # keep alkane-like peaks, sorted by RT
    alk = []
    for p in sorted(peaks, key=lambda x: x['apex_rt']):
        sc = data['scan_list'][p['apex_idx']]
        mz = np.asarray(sc['mz'], float); it = np.asarray(sc['intensity'], float)
        keep = mz >= 40                      # drop air for the alkane test
        mz, it = mz[keep], it[keep]
        if len(it) == 0 or it.max() <= 0:
            continue
        if not is_alkane(mz, it):
            continue
        # molecular ion candidate: highest m/z with >=2% rel intensity matching 14n+2
        rel = it / it.max() * 100
        mplus = None
        for j in np.argsort(mz)[::-1]:
            m = int(round(mz[j]))
            if rel[j] >= 1.5 and (m - 2) % 14 == 0 and 140 <= m <= 300:
                mplus = m; break
        n_from_m = (mplus - 2) // 14 if mplus else None
        alk.append({'rt': round(p['apex_rt'], 3), 'apex': p['apex_intensity'],
                    'mplus': mplus, 'n': n_from_m})

    # keep the strongest ~20 (the real alkane peaks dominate)
    alk = sorted(alk, key=lambda x: -x['apex'])[:25]
    alk = sorted(alk, key=lambda x: x['rt'])
    print(f"\n{'RT':>7} {'apexInt':>12} {'M+':>5} {'Cn(fromM+)':>10}")
    for a in alk:
        print(f"{a['rt']:7.2f} {a['apex']:12.0f} {str(a['mplus'] or '-'):>5} {str(a['n'] or '-'):>10}")


if __name__ == "__main__":
    main()
