#!/usr/bin/env python3
"""
GC-MS Auto-Processing Pipeline v2.0
====================================
10-Layer workflow per GCMS_Complete_Workflow_v2.md.
Key speed improvements: SNIP (vectorized), sparse representation,
centroid data, spectrum enhancement.
"""
import sys, json, argparse, time
from pathlib import Path

import numpy as np
import pandas as pd

from config import *
from step1_parse import convert_raw_to_mzml, load_mzml_to_matrix
from step2_preprocess import (preprocess_signal, estimate_column_bleed,
                               subtract_column_bleed)
from step3_peak_detect import detect_chromatographic_peaks
from step4_deconvolution import build_intensity_matrix, deconvolve_all_peaks
from step4_enhance import enhance_all_peaks, enhance_spectrum
from step5_integration import integrate_all_peaks
from step6_rt_align import align_retention_times
from step7_library_search import (SpectralLibrary, search_library, calc_ri_from_rt)
from step8_quantification import QuantificationEngine
from step9_export import compile_results, export_to_excel
from step10_qc import batch_qc_check


def process_single_sample(mzml_path: str, lib: SpectralLibrary,
                           config: dict = None,
                           output_dir: Path = None) -> dict:
    """Run v2 pipeline layers 0-8 on one sample."""
    cfg = dict(config or {})
    sample_name = Path(mzml_path).stem
    print(f"\n{'='*55}")
    print(f"  {sample_name}")
    print(f"{'='*55}")

    # ---- Layer 1: Load mzML ----
    t0 = time.time()
    print("  [L1] Loading mzML...")
    data = load_mzml_to_matrix(mzml_path)
    print(f"       {data['n_scans']:,} scans, "
          f"RT {data['rt'][0]:.1f}-{data['rt'][-1]:.1f} min")

    # ---- Build intensity matrix (sparse-friendly) ----
    intensity_matrix, mz_bins = build_intensity_matrix(
        data['scan_list'], data['rt'],
        mz_min=cfg.get('mz_min', MZ_MIN),
        mz_max=cfg.get('mz_max', MZ_MAX),
        mz_step=cfg.get('mz_step', MZ_STEP)
    )

    # ---- Layer 0: Column bleed subtraction ----
    print("  [L0] Column bleed subtraction...")
    bleed = estimate_column_bleed(intensity_matrix)
    intensity_matrix = subtract_column_bleed(intensity_matrix, bleed, scale=0.5)
    tic_raw = intensity_matrix.sum(axis=1)

    # ---- Layer 2: Preprocessing (S-G + SNIP) ----
    print("  [L2] S-G smooth + SNIP baseline...")
    pre = preprocess_signal(
        tic_raw,
        window_length=cfg.get('smooth_window', SMOOTH_WINDOW),
        max_half_window=cfg.get('snip_half_window', 100)
    )
    signal = pre['corrected']

    # ---- Layer 3: Peak detection (ICIS) ----
    print("  [L3] ICIS peak detection...")
    peaks = detect_chromatographic_peaks(
        signal, data['rt'],
        min_sn=cfg.get('min_sn', MIN_SN),
        min_peak_width_scans=cfg.get('min_peak_width', MIN_PEAK_WIDTH_SCANS),
        solvent_delay_min=4.0
    )
    print(f"       {len(peaks)} peaks detected")

    if not peaks:
        return {'sample_name': sample_name, 'dataframe': None, 'error': 'No peaks'}

    # ---- Layer 4: Spectrum Enhancement (SKIPPED — degrades NIST match quality) ----
    print("  [L4] Spectrum enhancement: OFF (binned apex spectrum performs better)")
    # peaks = enhance_all_peaks(intensity_matrix, peaks)  # DISABLED

    # ---- Layer 5: Deconvolution (OFF by default, use --deconv to enable) ----
    deconv_enabled = cfg.get('deconv_enabled', DECONV_ENABLED)
    if deconv_enabled:
        print("  [L5] Deconvolution (--deconv enabled)...")
        # Valley ratio check: only deconvolve peaks with deep valley (co-elution)
        deconv_peaks = []
        for i, p in enumerate(peaks):
            if i == 0 or i == len(peaks)-1:
                deconv_peaks.append(p)
                continue
            # Check valley depth between this peak and neighbors
            valley_left = signal[p['start_idx']-1] if p['start_idx'] > 0 else signal[p['start_idx']]
            valley_right = signal[p['end_idx']+1] if p['end_idx'] < len(signal)-1 else signal[p['end_idx']]
            valley = min(valley_left, valley_right)
            valley_ratio = valley / max(p['apex_intensity'], 1)
            if valley_ratio > cfg.get('deconv_valley_ratio', DECONV_VALLEY_RATIO):
                deconv_peaks.append(p)  # Deep valley = possible co-elution
        deconv = deconvolve_all_peaks(
            intensity_matrix, mz_bins, deconv_peaks,
            correlation_threshold=cfg.get('corr_threshold', CORRELATION_THRESHOLD)
        )
        # Map back: peaks NOT in deconv_peaks get empty deconv result
        deconv_full = []
        di = 0
        for p in peaks:
            if any(dp['apex_idx'] == p['apex_idx'] for dp in deconv_peaks):
                deconv_full.append(deconv[di] if di < len(deconv) else [])
                di += 1
            else:
                deconv_full.append([])
        deconv = deconv_full
    else:
        print("  [L5] Deconvolution: OFF (use --deconv to enable)")
        deconv = [[] for _ in peaks]

    # ---- Compute RI for all peaks first (needed for both search modes) ----
    from step7_library_search import calc_ri_from_rt
    for peak in peaks:
        peak['ri_measured'] = calc_ri_from_rt(peak['apex_rt'], ALKANE_RTS_DB_WAX)

    # ---- Layer 6: Library search ----
    use_nist = cfg.get('use_nist', False)
    if use_nist:
        print("  [L6] NIST MS Search engine + RI re-ranking...")
        from nist_engine import NISTSearchEngine
        nist_engine = NISTSearchEngine()
        # Use RAW spectrum (proven 29% coverage vs 24% binned)
        id_results = nist_engine.search_peaks_raw(
            peaks, data['scan_list'], top_n=5, verbose=True
        )
        n_confident = sum(1 for m in id_results if m and m[0].get('rmf', 0) >= 800)
    else:
        print("  [L6] Library search + RI boost...")
        min_rmf = cfg.get('min_rmf', MIN_RMF)
        ri_tol = cfg.get('ri_tolerance', RI_TOLERANCE)

    if not use_nist:
        id_results = []
        n_confident = 0
        for i, peak in enumerate(peaks):
            if deconv_enabled and deconv[i]:
                query = deconv[i][0]['pure_spectrum']
            elif peak.get('enhanced_full') is not None:
                query = peak['enhanced_full']
            else:
                query = intensity_matrix[peak['apex_idx'], :]
            peak['ri_measured'] = ri_meas

            matches = search_library(
                query, lib, mz_bins, top_n=10, max_candidates=500,
                min_rmf=min_rmf, ri_tolerance=ri_tol,
                query_ri=ri_meas, ri_weight=RI_WEIGHT,
                coelution_penalty=True
            )
            id_results.append(matches)
            peak['ri_lit'] = matches[0].get('ri') if matches else None

            if matches and matches[0].get('rmf', 0) >= 800:
                n_confident += 1

    print(f"       {n_confident}/{len(peaks)} peaks matched (RMF>=800)")

    # ---- Layer 7: Integration ----
    print("  [L7] ICIS integration...")
    integrated = integrate_all_peaks(signal, data['rt'], peaks,
                                      use_peak_baseline=USE_PEAK_BASELINE)

    # ---- Export MSP for NIST MS Search ----
    msp_path = output_dir / f"{sample_name}_spectra.msp"
    from nist_bridge import export_peaks_to_msp
    export_peaks_to_msp(integrated, intensity_matrix, mz_bins, data['rt'],
                         str(msp_path), sample_name)
    nist_msp_path = str(msp_path)

    # ---- Layer 8: Quantification ----
    qe = QuantificationEngine()
    all_areas = [p['area'] for p in integrated]
    quant_results = [
        qe.area_normalization(p['area'], all_areas)
        if cfg.get('quant_method', QUANT_METHOD) == 'normalization'
        else {'method': cfg.get('quant_method', QUANT_METHOD)}
        for p in integrated
    ]

    # ---- Compile ----
    result_df = compile_results(sample_name, integrated, id_results, quant_results)
    elapsed = time.time() - t0
    print(f"       Done in {elapsed:.0f}s")

    return {
        'sample_name': sample_name,
        'dataframe': result_df,
        'integrated_peaks': integrated,
        'data': data,
        'nist_msp': nist_msp_path
    }


def run_gcms_pipeline(raw_files: list = None, mzml_files: list = None,
                       output_dir: str = None, library_path: str = None,
                       config: dict = None) -> dict:
    """Main v2 pipeline entry point."""
    if output_dir is None:
        output_dir = OUTPUT_DIR
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    if library_path is None:
        library_path = LIBRARY_PATH
    cfg = dict(config or {})

    print("=" * 55)
    print("  GC-MS Auto-Processing Pipeline v2.0")
    print("  SNIP baseline | ICIS peaks | Enhanced spectra | Triple ID")
    print("=" * 55)

    # ---- Load library ----
    print(f"\n[LIB] Loading: {Path(library_path).name}")
    mz_bins_load = np.arange(MZ_MIN, MZ_MAX + MZ_STEP, MZ_STEP)
    lib = SpectralLibrary()
    lib.load_json(str(library_path), mz_bins_load)

    # ---- Step 1: RAW conversion ----
    if mzml_files is None:
        mzml_files = []
    if raw_files:
        print(f"\n[L0] Converting {len(raw_files)} RAW files...")
        for raw_f in raw_files:
            mzml_f = convert_raw_to_mzml(raw_f, str(output_dir / "mzml"))
            mzml_files.append(mzml_f)
            print(f"     {Path(raw_f).name} → {Path(mzml_f).name}")

    # ---- Process each sample ----
    all_results = {}
    for mzml_f in mzml_files:
        result = process_single_sample(mzml_f, lib, cfg, output_dir)
        if result['dataframe'] is not None:
            all_results[result['sample_name']] = result['dataframe']

    # ---- Layer 9: Cross-sample alignment ----
    if len(all_results) > 1:
        print(f"\n[L9] RT alignment across {len(all_results)} samples...")
        ref_name = list(all_results.keys())[0]
        ref_peaks = [{'apex_rt': r['RT_min']}
                     for _, r in all_results[ref_name].iterrows()
                     if r['Compound_Name'] != 'Unknown']
        for sn in all_results:
            if sn == ref_name:
                continue
            sp = [{'apex_rt': r['RT_min']}
                  for _, r in all_results[sn].iterrows()
                  if r['Compound_Name'] != 'Unknown']
            aligned = align_retention_times(sp, ref_peaks,
                                             rt_tolerance=cfg.get('rt_align_tol',
                                                                  RT_ALIGN_TOLERANCE))
            rt_map = {a['apex_rt_original']: a['apex_rt']
                      for a in aligned if 'apex_rt_original' in a}
            for idx, row in all_results[sn].iterrows():
                if row['RT_min'] in rt_map:
                    all_results[sn].at[idx, 'RT_min'] = round(rt_map[row['RT_min']], 3)
        print("       Done")

    # ---- Layer 10: QC audit ----
    print(f"\n[L10] QC audit...")
    qc = batch_qc_check(all_results,
                        qc_sample_pattern=cfg.get('qc_pattern', QC_SAMPLE_PATTERN))
    print(f"       Status: {qc['status']}")
    for w in qc.get('warnings', []):
        print(f"       [WARN] {w}")

    # ---- Export ----
    if all_results:
        combined = pd.concat(all_results.values(), ignore_index=True)
        # Use fixed filename based on first sample name (overwrites previous runs)
        first_sample = list(all_results.keys())[0] if all_results else "results"
        xlsx = str(output_dir / f"{first_sample}.xlsx")
        export_to_excel(combined, xlsx)
    else:
        xlsx = None

    return {
        'results': all_results, 'qc_report': qc,
        'output_file': xlsx,
        'total': len(mzml_files), 'processed': len(all_results)
    }


# ---- CLI ----
if __name__ == "__main__":
    ap = argparse.ArgumentParser(description='GC-MS Pipeline v2.0')
    ap.add_argument('--raw-files', nargs='+', help='Input .RAW files')
    ap.add_argument('--mzml-files', nargs='+', help='Input .mzML files')
    ap.add_argument('--output', '-o', default=None, help='Output directory')
    ap.add_argument('--library', '-l', default=None, help='Library JSON')
    ap.add_argument('--min-sn', type=float, default=None)
    ap.add_argument('--min-rmf', type=int, default=None)
    ap.add_argument('--deconv', action='store_true', help='Enable deconvolution')
    ap.add_argument('--nist', action='store_true', help='Use NIST MS Search engine (official)')
    ap.add_argument('--config', default=None, help='JSON config override')
    args = ap.parse_args()

    if not args.raw_files and not args.mzml_files:
        ap.error("Specify --raw-files or --mzml-files")

    config = {}
    if args.config:
        with open(args.config) as f:
            config = json.load(f)
    if args.min_sn: config['min_sn'] = args.min_sn
    if args.min_rmf: config['min_rmf'] = args.min_rmf
    if args.deconv: config['deconv_enabled'] = True
    if args.nist: config['use_nist'] = True

    result = run_gcms_pipeline(
        raw_files=args.raw_files, mzml_files=args.mzml_files,
        output_dir=args.output, library_path=args.library,
        config=config if config else None
    )

    print(f"\n{'='*55}")
    print(f"  Pipeline complete: {result['processed']}/{result['total']} samples")
    if result['output_file']:
        print(f"  Output: {result['output_file']}")
    print(f"{'='*55}")
