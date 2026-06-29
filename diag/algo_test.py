#!/usr/bin/env python3
"""
diag/algo_test.py — is the discrepancy an ALGORITHM problem or a SPECTRUM problem?
================================================================================
Two controlled tests:

A) SELF-SEARCH FAITHFULNESS: take spectra straight from the NIST library and
   search them back through the pipeline's exact call. A correct algorithm MUST
   return the same compound at RMF ~999. If it doesn't, the call is broken.

B) NOMINAL-MASS BINNING: NIST EI search is nominal-mass. The pipeline currently
   sends raw float m/z from the mzML. Re-score the 49 with proper integer
   nominal-mass binning (round m/z, sum intensities per unit mass). If coverage
   improves, the float-m/z handling is a real algorithm defect.
"""
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "eval"))

from pyms.Spectrum import MassSpectrum
from config import (MZ_MIN, MZ_MAX, MZ_STEP, SMOOTH_WINDOW, MIN_SN,
                    MIN_PEAK_WIDTH_SCANS, ALKANE_RTS_DB_WAX)
from step1_parse import load_mzml_to_matrix
from step2_preprocess import estimate_column_bleed, subtract_column_bleed, preprocess_signal
from step3_peak_detect import detect_chromatographic_peaks
from step4_deconvolution import build_intensity_matrix
from step7_library_search import calc_ri_from_rt
from nist_engine import NISTSearchEngine
from score import load_ground_truth, score, DEFAULT_GROUND_TRUTH
from bm_probe import MZML


def nominal_bin(mz, inten):
    """Collapse to unit (nominal) mass: round m/z, sum intensities per integer."""
    mz = np.asarray(mz, float); inten = np.asarray(inten, float)
    d = {}
    for m, i in zip(mz, inten):
        k = int(round(m))
        d[k] = d.get(k, 0.0) + i
    ks = sorted(d)
    return np.array(ks, float), np.array([d[k] for k in ks], float)


def main():
    eng = NISTSearchEngine()
    data = load_mzml_to_matrix(MZML)

    # ---------- Test A: self-search faithfulness ----------
    print("=" * 64)
    print("  TEST A — search NIST library spectra against the library")
    print("=" * 64)
    sc = data['scan_list'][2000]
    hits = eng.engine.full_search_with_ref_data(
        MassSpectrum(*nominal_bin(sc['mz'], sc['intensity'])))
    checked = 0
    for sr, ref in hits[:6]:
        ms = ref.mass_spec
        back = eng.engine.full_search_with_ref_data(
            MassSpectrum(list(ms.mass_list), list(ms.intensity_list)))
        bsr, bref = back[0]
        ok = bref.name == ref.name
        print(f"  {ref.name[:42]:<42} self-RMF={bsr.reverse_match_factor:>3} "
              f"{'OK' if ok else 'MISMATCH'}")
        checked += 1
    print("  -> RMF~999 + OK means the NIST call is faithful (algorithm correct)")

    # ---------- Test B: nominal binning on the full 49 ----------
    print("\n" + "=" * 64)
    print("  TEST B — float m/z (current) vs nominal-mass binning, on the 49")
    print("=" * 64)
    matrix, mz_bins = build_intensity_matrix(data['scan_list'], data['rt'],
                                             MZ_MIN, MZ_MAX, MZ_STEP)
    matrix = subtract_column_bleed(matrix, estimate_column_bleed(matrix), scale=0.5)
    pre = preprocess_signal(matrix.sum(axis=1), window_length=SMOOTH_WINDOW, max_half_window=100)
    peaks = detect_chromatographic_peaks(pre['corrected'], data['rt'], min_sn=MIN_SN,
                                         min_peak_width_scans=MIN_PEAK_WIDTH_SCANS,
                                         solvent_delay_min=4.0)
    for p in peaks:
        p['ri_measured'] = calc_ri_from_rt(p['apex_rt'], ALKANE_RTS_DB_WAX)

    def build(use_nominal):
        cand = []
        for p in peaks:
            sc = data['scan_list'][p['apex_idx']]
            mz, it = (nominal_bin(sc['mz'], sc['intensity']) if use_nominal
                      else (sc['mz'], sc['intensity']))
            res = eng.search_raw_spectrum(np.asarray(mz, float), np.asarray(it, float),
                                          top_n=5, measured_ri=p['ri_measured'])
            cand.append({'rt': round(p['apex_rt'], 3),
                         'candidates': [{'name': m['name'], 'rmf': m['rmf']} for m in res]})
        return cand

    user = load_ground_truth(DEFAULT_GROUND_TRUTH)
    for label, nominal in (("float m/z (current)", False), ("nominal-mass binned", True)):
        res = score(user, build(nominal), top_n=5)
        t1 = sum(1 for r in res if r['top1']); t5 = sum(1 for r in res if r['top5'])
        print(f"  {label:<22}  top-1={t1}/49  top-5={t5}/49")


if __name__ == "__main__":
    main()
