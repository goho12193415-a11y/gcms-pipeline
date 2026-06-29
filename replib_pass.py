#!/usr/bin/env python3
"""
replib_pass.py — standalone NIST replicate-library (replib) search.
===================================================================
Runs in its OWN process because pyms_nist_search wraps a stateful DLL with a
single global active library: a replib Engine and a mainlib Engine cannot
coexist in one process (the second hijacks the first). The pipeline calls this
as a subprocess and merges the result offline (reserve strategy).

Rebuilds the exact same peaks as the pipeline (deterministic L0-L3) so the
output aligns to the pipeline's peaks by index.

    python replib_pass.py <input.mzML> <output.json>
"""
import sys, json
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import pyms_nist_search
from config import (MZ_MIN, MZ_MAX, MZ_STEP, SMOOTH_WINDOW, MIN_SN,
                    MIN_PEAK_WIDTH_SCANS, ALKANE_RTS_DB_WAX, TOP_N_CANDIDATES,
                    NIST_REPLIB_CANDIDATES)
from step1_parse import load_sample
from step2_preprocess import (preprocess_signal, estimate_column_bleed,
                              subtract_column_bleed)
from step3_peak_detect import detect_chromatographic_peaks
from step4_deconvolution import build_intensity_matrix
from step7_library_search import calc_ri_from_rt
from nist_engine import NISTSearchEngine


def main():
    mzml_path, out_path = sys.argv[1], sys.argv[2]
    replib = next((rp for rp in NIST_REPLIB_CANDIDATES if Path(rp).exists()), None)
    if replib is None:
        json.dump([], open(out_path, 'w', encoding='utf-8'))
        print("  [replib_pass] no replib found")
        return

    data = load_sample(mzml_path)
    matrix, mz_bins = build_intensity_matrix(data['scan_list'], data['rt'],
                                             MZ_MIN, MZ_MAX, MZ_STEP)
    matrix = subtract_column_bleed(matrix, estimate_column_bleed(matrix), scale=0.5)
    pre = preprocess_signal(matrix.sum(axis=1), window_length=SMOOTH_WINDOW,
                            max_half_window=100)
    peaks = detect_chromatographic_peaks(
        pre['corrected'], data['rt'], min_sn=MIN_SN,
        min_peak_width_scans=MIN_PEAK_WIDTH_SCANS, solvent_delay_min=4.0)
    for p in peaks:
        p['ri_measured'] = calc_ri_from_rt(p['apex_rt'], ALKANE_RTS_DB_WAX)

    eng = NISTSearchEngine(lib_path=replib,
                           lib_type=pyms_nist_search.NISTMS_REP_LIB)
    matches = eng.search_peaks_raw(peaks, data['scan_list'],
                                   top_n=TOP_N_CANDIDATES, verbose=False)

    dump = [{'peak_no': i + 1, 'rt': round(peaks[i]['apex_rt'], 3),
             'candidates': matches[i] or []}
            for i in range(len(peaks))]
    json.dump(dump, open(out_path, 'w', encoding='utf-8'), ensure_ascii=False)
    print(f"  [replib_pass] {len(peaks)} peaks searched against replib")


if __name__ == "__main__":
    main()
