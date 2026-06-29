#!/usr/bin/env python3
"""
diag/a_ri_inject.py — feasibility of RI-primary injection for A-class.
=====================================================================
For each A-class branched alkane:
  1. measured RI (from peak RT)
  2. is the EXACT compound in the predicted-RI DB (ri_nist_full.name_to_ri)?
     with what predicted RI, and |pred - measured|?
  3. among ALL branched alkanes in the predicted DB within +/-tol of measured
     RI, what rank would the correct one get (RI-primary candidate list)?

This tells us whether injecting RI-matched branched alkanes would surface the
right answer, and how bad the isomer ambiguity is.
"""
import sys, json, re
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "eval"))
from config import ALKANE_RTS_DB_WAX
from step7_library_search import calc_ri_from_rt
from score import normalize_alkane, names_match
from a_rank import A_TARGETS
from bm_probe import find_peak
from step1_parse import load_mzml_to_matrix
from step2_preprocess import estimate_column_bleed, subtract_column_bleed, preprocess_signal
from step3_peak_detect import detect_chromatographic_peaks
from step4_deconvolution import build_intensity_matrix
from config import MZ_MIN, MZ_MAX, MZ_STEP, SMOOTH_WINDOW, MIN_SN, MIN_PEAK_WIDTH_SCANS
from bm_probe import MZML

LIB = Path(__file__).resolve().parent.parent / "library" / "ri_nist_full.json"

ALKANE_RE = re.compile(r'(methyl|trimethyl|tetramethyl|dimethyl).*(ane)\b|(ane),.*methyl', re.I)


def is_branched_alkane(name):
    n = name.lower()
    return ('methyl' in n) and ('ane' in n) and 'ene' not in n and 'ol' not in n \
        and 'one' not in n and 'acid' not in n and 'ester' not in n


def main():
    d = json.load(open(LIB, encoding='utf-8'))
    name_to_ri = d['name_to_ri']
    # normalized index of branched alkanes in predicted DB
    branched = {}   # canon -> (name, ri)
    for nm, ri in name_to_ri.items():
        if is_branched_alkane(nm):
            branched[normalize_alkane(nm)] = (nm, ri)
    print(f"predicted-RI DB: {len(name_to_ri)} names, "
          f"{len(branched)} branched alkanes\n")

    # rebuild peaks for measured RI
    data = load_mzml_to_matrix(MZML)
    matrix, mz_bins = build_intensity_matrix(data['scan_list'], data['rt'],
                                             MZ_MIN, MZ_MAX, MZ_STEP)
    matrix = subtract_column_bleed(matrix, estimate_column_bleed(matrix), scale=0.5)
    pre = preprocess_signal(matrix.sum(axis=1), window_length=SMOOTH_WINDOW, max_half_window=100)
    peaks = detect_chromatographic_peaks(pre['corrected'], data['rt'], min_sn=MIN_SN,
                                         min_peak_width_scans=MIN_PEAK_WIDTH_SCANS,
                                         solvent_delay_min=4.0)

    print("=" * 92)
    print(f"  {'compound':<32} {'measRI':>7} {'predRI':>7} {'|d|':>5} {'inDB':>5} "
          f"{'rank/inWin':>11}")
    print("=" * 92)
    for tol in (30,):
        for rt, name in A_TARGETS:
            pk = find_peak(peaks, rt)
            mri = calc_ri_from_rt(pk['apex_rt'], ALKANE_RTS_DB_WAX) if pk else None
            canon = normalize_alkane(name)
            hit = branched.get(canon)
            pred = hit[1] if hit else None
            dd = abs(pred - mri) if (pred and mri) else None
            # candidates within window of measured RI
            win = sorted(((abs(ri - mri), nm, ri) for c, (nm, ri) in branched.items()
                          if mri and abs(ri - mri) <= tol), key=lambda x: x[0])
            rank = None
            for k, (_, nm, _) in enumerate(win):
                if names_match(name, nm):
                    rank = k + 1
                    break
            rankstr = f"{rank}/{len(win)}" if rank else f"-/{len(win)}"
            print(f"  {name:<32} {mri or 0:7.0f} {pred or 0:7.0f} "
                  f"{(dd if dd is not None else -1):5.0f} {('Y' if hit else '-'):>5} {rankstr:>11}")
    print("=" * 92)
    print("  inDB=exact compound has a predicted RI;  rank/inWin = its rank among")
    print("  branched alkanes within +/-30 RI of measured (isomer ambiguity size)")


if __name__ == "__main__":
    main()
