#!/usr/bin/env python3
"""
diag/validate_qgd.py — validate the pipeline on a Shimadzu .qgd vs its
GCMSsolution result table. Cross-instrument generalization test.

    python diag/validate_qgd.py <sample.qgd> <result.xlsx>
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "eval"))

from config import MZ_MIN, MZ_MAX, MZ_STEP, SMOOTH_WINDOW, MIN_SN, MIN_PEAK_WIDTH_SCANS
from step1_parse import load_sample
from step2_preprocess import estimate_column_bleed, subtract_column_bleed, preprocess_signal
from step3_peak_detect import detect_chromatographic_peaks
from step4_deconvolution import build_intensity_matrix
from nist_engine import NISTSearchEngine
from score import names_match


def load_gt(xlsx):
    df = pd.read_excel(xlsx, sheet_name=0, header=8)
    df = df[pd.to_numeric(df['Peak#'], errors='coerce').notna()]
    out = []
    for _, r in df.iterrows():
        rt = pd.to_numeric(r['Ret.Time'], errors='coerce')
        nm = str(r['Name']).strip()
        if pd.notna(rt) and nm and nm != 'nan':
            out.append({'rt': float(rt), 'name': nm, 'si': r.get('SI', '')})
    return out


def main():
    qgd, xlsx = sys.argv[1], sys.argv[2]
    data = load_sample(qgd)
    print(f"loaded {data['n_scans']} scans, RT {data['rt'][0]:.1f}-{data['rt'][-1]:.1f}, "
          f"m/z {data['mz_range'][0]:.0f}-{data['mz_range'][1]:.0f}")
    matrix, mz_bins = build_intensity_matrix(data['scan_list'], data['rt'], MZ_MIN, MZ_MAX, MZ_STEP)
    matrix = subtract_column_bleed(matrix, estimate_column_bleed(matrix), scale=0.5)
    pre = preprocess_signal(matrix.sum(axis=1), window_length=SMOOTH_WINDOW, max_half_window=100)
    peaks = detect_chromatographic_peaks(pre['corrected'], data['rt'], min_sn=MIN_SN,
                                         min_peak_width_scans=MIN_PEAK_WIDTH_SCANS,
                                         solvent_delay_min=2.0)
    print(f"{len(peaks)} peaks detected")
    eng = NISTSearchEngine()
    matches = eng.search_peaks_raw(peaks, data['scan_list'], top_n=8, verbose=False)
    cand = [{'rt': round(peaks[i]['apex_rt'], 3),
             'cands': [m.get('name', '') for m in (matches[i] or [])]} for i in range(len(peaks))]

    gt = load_gt(xlsx)
    def peak_at(rt, tol=0.1):
        best, bd = None, tol
        for c in cand:
            d = abs(c['rt'] - rt)
            if d < bd:
                bd, best = d, c
        return best

    n = len(gt); t1 = t5 = t8 = nopk = 0
    misses = []; hits = []
    for g in gt:
        pk = peak_at(g['rt'])
        if pk is None:
            nopk += 1; misses.append((g, None)); continue
        cs = pk['cands']
        h1 = bool(cs) and names_match(g['name'], cs[0])
        h5 = any(names_match(g['name'], c) for c in cs[:5])
        h8 = any(names_match(g['name'], c) for c in cs[:8])
        if h1: t1 += 1
        if h5: t5 += 1
        if h8: t8 += 1; hits.append((g, [c for c in cs if names_match(g['name'], c)][:1]))
        else: misses.append((g, pk))
    print("=" * 66)
    print(f"  vs GCMSsolution result ({n} peaks), mainlib top-N, RT tol 0.1:")
    print(f"    no pipeline peak within tol: {nopk}")
    print(f"    top-1: {t1}/{n} = {t1*100//n}%")
    print(f"    top-5: {t5}/{n} = {t5*100//n}%")
    print(f"    top-8: {t8}/{n} = {t8*100//n}%")
    print(f"  HITS (top-8, {len(hits)}):")
    for g, _ in sorted(hits, key=lambda z: z[0]['rt']):
        print(f"    RT={g['rt']:5.2f} SI={g['si']} {g['name'][:38]}")
    print("=" * 66)
    print("  misses (GCMSsolution name not in pipeline top-8):")
    for g, pk in sorted(misses, key=lambda z: z[0]['rt']):
        top = pk['cands'][0] if (pk and pk['cands']) else '[no peak]'
        print(f"    RT={g['rt']:5.2f} SI={g['si']} {g['name'][:34]:<34} -> {top[:30]}")


if __name__ == "__main__":
    main()
