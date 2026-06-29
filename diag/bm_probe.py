#!/usr/bin/env python3
"""
diag/bm_probe.py — Bm-class root-cause diagnosis.
=================================================
The archive proposes 3 mechanisms for the Bm-class failures
(library HAS the compound, but pipeline top-5 misses it):
  M1  peak boundary wrong  -> neighbor-ion contamination
  M2  apex scan wrong      -> TIC max is not the purest scan
  M3  mzML conversion drop -> low-abundance characteristic ions lost

This probe rebuilds the pipeline state and, for each Bm peak, sends THREE
spectrum variants to the SAME NIST engine and reports whether the correct
answer comes back in the top-5:
  A  apex single raw scan         (current pipeline method)
  B  apex +/-3 scans, averaged    (tests M2 / recovers low-abundance ions)
  C  co-elution-aware component   (tests M1 / removes neighbor ions)

It also prints the RT gap to the nearest neighbor peaks (M1 evidence).

Usage:  python diag/bm_probe.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from config import (MZ_MIN, MZ_MAX, MZ_STEP, SMOOTH_WINDOW, MIN_SN,
                    MIN_PEAK_WIDTH_SCANS, ALKANE_RTS_DB_WAX)
from step1_parse import load_mzml_to_matrix
from step2_preprocess import (preprocess_signal, estimate_column_bleed,
                              subtract_column_bleed)
from step3_peak_detect import detect_chromatographic_peaks
from step4_deconvolution import build_intensity_matrix
from step7_library_search import calc_ri_from_rt
from nist_engine import NISTSearchEngine

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "eval"))
from score import names_match

MZML = r"C:\gcms_pipeline\output\mzml\xt6.26.mzML"

# Archive's 8 Bm-class compounds (RT, expected name)
BM_TARGETS = [
    (12.82, "(E)-2-Heptenal"),
    (13.07, "2,3-Octanedione"),
    (19.39, "2,4-Dimethyl-3-ethylpyrrole"),
    (20.76, "(E)-2-Nonenal"),
    (22.98, "2,6,10-Trimethylpentadecane"),
    (23.53, "beta-Cyclocitral"),
    (26.90, "2,6-Di-tert-butyl-4-methylphenol"),  # BHT
    (37.22, "alpha-Isomethylionone"),
]


def find_peak(peaks, rt, tol=0.3):
    best, bestd = None, tol
    for p in peaks:
        d = abs(p['apex_rt'] - rt)
        if d < bestd:
            bestd, best = d, p
    return best


def spec_apex_raw(scan_list, peak):
    """Strategy A: single apex raw scan (current method)."""
    s = scan_list[peak['apex_idx']]
    return np.asarray(s['mz'], float), np.asarray(s['intensity'], float)


def spec_apex_avg(matrix, mz_bins, peak, half=3):
    """Strategy B: average binned spectra over apex +/- half scans."""
    a = peak['apex_idx']
    lo, hi = max(0, a - half), min(matrix.shape[0] - 1, a + half)
    spec = matrix[lo:hi + 1, :].mean(axis=0)
    mask = spec > 0
    return mz_bins[mask], spec[mask]


def spec_component(matrix, mz_bins, peak, corr_thresh=0.85):
    """Strategy C: co-elution-aware component spectrum.

    Keep only m/z bins whose in-peak EIC correlates with the peak's
    chromatographic profile; intensity taken at the apex. This removes
    ions that belong to co-eluting neighbors (M1)."""
    s, e, a = peak['start_idx'], peak['end_idx'], peak['apex_idx']
    s = max(0, s); e = min(matrix.shape[0] - 1, e)
    if e - s < 3:
        s, e = max(0, a - 4), min(matrix.shape[0] - 1, a + 4)
    region = matrix[s:e + 1, :]                 # (scans, bins)
    ref = region.sum(axis=1)                    # TIC profile in region
    ref = ref - ref.mean()
    ref_norm = np.linalg.norm(ref)
    if ref_norm == 0:
        return spec_apex_avg(matrix, mz_bins, peak)
    apex_local = a - s
    keep_mz, keep_int = [], []
    for j in range(region.shape[1]):
        col = region[:, j]
        if col.max() <= 0:
            continue
        cj = col - col.mean()
        nj = np.linalg.norm(cj)
        if nj == 0:
            continue
        corr = float(np.dot(cj, ref) / (nj * ref_norm))
        if corr >= corr_thresh:
            keep_mz.append(mz_bins[j])
            keep_int.append(region[apex_local, j])
    if not keep_mz:
        return spec_apex_avg(matrix, mz_bins, peak)
    return np.asarray(keep_mz, float), np.asarray(keep_int, float)


def in_top5(engine, mz, inten, expected, ri):
    res = engine.search_raw_spectrum(mz, inten, top_n=5, measured_ri=ri)
    hit_rank = None
    for k, r in enumerate(res):
        if names_match(expected, r['name']):
            hit_rank = k + 1
            break
    top1 = res[0]['name'] if res else '-'
    top1_rmf = res[0]['rmf'] if res else 0
    return hit_rank, top1, top1_rmf


def main():
    print("Loading mzML + rebuilding pipeline state...")
    data = load_mzml_to_matrix(MZML)
    matrix, mz_bins = build_intensity_matrix(
        data['scan_list'], data['rt'], MZ_MIN, MZ_MAX, MZ_STEP)
    bleed = estimate_column_bleed(matrix)
    matrix = subtract_column_bleed(matrix, bleed, scale=0.5)
    tic = matrix.sum(axis=1)
    pre = preprocess_signal(tic, window_length=SMOOTH_WINDOW, max_half_window=100)
    peaks = detect_chromatographic_peaks(
        pre['corrected'], data['rt'], min_sn=MIN_SN,
        min_peak_width_scans=MIN_PEAK_WIDTH_SCANS, solvent_delay_min=4.0)
    print(f"  {len(peaks)} peaks detected\n")

    engine = NISTSearchEngine()
    all_rts = sorted(p['apex_rt'] for p in peaks)

    print("=" * 78)
    print(f"  {'RT':>6} {'compound':<28} {'A apex':>7} {'B avg':>7} {'C comp':>7}  nbr-dRT")
    print("=" * 78)
    summary = {'A': 0, 'B': 0, 'C': 0}
    for rt, name in BM_TARGETS:
        pk = find_peak(peaks, rt)
        if pk is None:
            print(f"  {rt:6.2f} {name:<28} [no peak]")
            continue
        ri = calc_ri_from_rt(pk['apex_rt'], ALKANE_RTS_DB_WAX)

        mzA, inA = spec_apex_raw(data['scan_list'], pk)
        mzB, inB = spec_apex_avg(matrix, mz_bins, pk)
        mzC, inC = spec_component(matrix, mz_bins, pk)

        rA, t1A, rmfA = in_top5(engine, mzA, inA, name, ri)
        rB, _, _ = in_top5(engine, mzB, inB, name, ri)
        rC, _, _ = in_top5(engine, mzC, inC, name, ri)

        # nearest-neighbor RT gap (M1 evidence)
        idx = all_rts.index(pk['apex_rt'])
        left = pk['apex_rt'] - all_rts[idx - 1] if idx > 0 else 9.9
        right = all_rts[idx + 1] - pk['apex_rt'] if idx < len(all_rts) - 1 else 9.9
        nbr = min(left, right)

        def fmt(r):
            return f"#{r}" if r else "-"
        if rA:
            summary['A'] += 1
        if rB:
            summary['B'] += 1
        if rC:
            summary['C'] += 1
        flag = " <co-elute" if nbr < 0.10 else ""
        print(f"  {rt:6.2f} {name:<28} {fmt(rA):>7} {fmt(rB):>7} {fmt(rC):>7}  {nbr:.3f}{flag}")
        print(f"         apex top-1: {t1A[:50]:<50} RMF={rmfA}")

    print("=" * 78)
    n = len(BM_TARGETS)
    print(f"  Recovered in top-5:  A(apex)={summary['A']}/{n}  "
          f"B(avg)={summary['B']}/{n}  C(component)={summary['C']}/{n}")
    print("  '#k' = correct answer found at rank k;  '-' = not in top-5")
    print("  nbr-dRT<0.10 flags possible co-elution (mechanism M1)")


if __name__ == "__main__":
    main()
