"""
ri_calib.py — derive an n-alkane RI ladder from a standard run.
===============================================================
Per-sample / per-instrument RI calibration: given an n-alkane standard
(.qgd / .raw / .mzML), detect the alkane peaks, assign carbon number from the
molecular ion (M+ = 14n+2), and return a clean monotonic RT->RI ladder
(RI = 100*n) usable by step7_library_search.calc_ri_from_rt.

This replaces the single hard-coded Thermo DB-WAX ladder so each sample can be
calibrated against the standard actually run on its instrument/column.
"""
import numpy as np

from config import (MZ_MIN, MZ_MAX, MZ_STEP, SMOOTH_WINDOW, MIN_SN,
                    MIN_PEAK_WIDTH_SCANS)
from step1_parse import load_sample
from step2_preprocess import (estimate_column_bleed, subtract_column_bleed,
                              preprocess_signal)
from step3_peak_detect import detect_chromatographic_peaks
from step4_deconvolution import build_intensity_matrix


def _is_alkane(mz, it):
    if len(mz) == 0:
        return False
    base = int(round(mz[np.argmax(it)]))
    mzr = np.round(mz).astype(int)
    series = sum(1 for s in (43, 57, 71, 85, 99) if s in mzr)
    return base in (43, 57, 71, 85) and series >= 4


def build_ri_ladder(standard_path, solvent_delay=1.5, verbose=True):
    """Return a clean monotonic [(rt, RI), ...] ladder from an alkane standard."""
    data = load_sample(standard_path)
    matrix, mz_bins = build_intensity_matrix(data['scan_list'], data['rt'],
                                             MZ_MIN, MZ_MAX, MZ_STEP)
    matrix = subtract_column_bleed(matrix, estimate_column_bleed(matrix), scale=0.5)
    pre = preprocess_signal(matrix.sum(axis=1), window_length=SMOOTH_WINDOW,
                            max_half_window=100)
    peaks = detect_chromatographic_peaks(pre['corrected'], data['rt'], min_sn=MIN_SN,
                                         min_peak_width_scans=MIN_PEAK_WIDTH_SCANS,
                                         solvent_delay_min=solvent_delay)
    # collect ALL alkane candidates (rt, carbon#) identified by M+
    cand = []
    for p in peaks:
        sc = data['scan_list'][p['apex_idx']]
        mz = np.asarray(sc['mz'], float); it = np.asarray(sc['intensity'], float)
        k = mz >= 40
        mz, it = mz[k], it[k]
        if not _is_alkane(mz, it):
            continue
        rel = it / it.max() * 100
        mplus = None
        for j in np.argsort(mz)[::-1]:
            m = int(round(mz[j]))
            if rel[j] >= 1.0 and (m - 2) % 14 == 0 and 120 <= m <= int(MZ_MAX):
                mplus = m; break
        if mplus is None:
            continue
        n = (mplus - 2) // 14
        if 8 <= n <= 40:
            cand.append((p['apex_rt'], n))

    # the true ladder is monotonic (RT up <=> carbon# up); take the longest
    # increasing subsequence of carbon number after sorting by RT — this drops
    # spurious M+ matches automatically.
    cand.sort(key=lambda x: x[0])
    seq = [c[1] for c in cand]
    L = len(seq)
    if L == 0:
        if verbose:
            print("  [RI-calib] no alkane ladder found")
        return []
    best = [1] * L; prev = [-1] * L
    for i in range(L):
        for j in range(i):
            if seq[j] < seq[i] and best[j] + 1 > best[i]:
                best[i] = best[j] + 1; prev[i] = j
    end = int(np.argmax(best))
    chain = []
    while end != -1:
        chain.append(end); end = prev[end]
    chain.reverse()
    ladder = [(round(cand[i][0], 3), 100 * cand[i][1]) for i in chain]
    if verbose:
        print(f"  [RI-calib] {len(ladder)} alkanes  "
              f"C{ladder[0][1] // 100}-C{ladder[-1][1] // 100}")
    return ladder


if __name__ == "__main__":
    import sys
    lad = build_ri_ladder(sys.argv[1])
    for rt, ri in lad:
        print(f"  RT {rt:7.3f}  RI {ri}")
