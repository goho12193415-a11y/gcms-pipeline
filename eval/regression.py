#!/usr/bin/env python3
"""
eval/regression.py — regression guard + verified anchors (+ optional accuracy).
==============================================================================
Runs the pipeline on the reference sample and checks three things:

  1. AGREEMENT vs the 49 manual top-1 picks (top-1/5/8).
     *** REGRESSION-ONLY METRIC, NOT ACCURACY ***
     The manual reference is itself fallible (it is just NIST top-1 picks and
     contains known errors), so this measures "did a change move the pipeline
     toward/away from the established manual result", NOT correctness. Floors
     exist only to catch regressions.

  2. VERIFIED ANCHORS must all stay found (top-8). These 15 are the dual-
     evidence (strong MS + RI-on-calibration-line) confident IDs. Losing any
     one is a hard FAIL. (Note: this is a must-not-regress floor, not an
     accuracy measure — the anchors were defined using the pipeline's own
     strong hits, so it is circular as an accuracy denominator.)

  3. ACCURACY (optional, --truth-set FILE): the ONLY real accuracy number.
     Provide a file of compounds known a-priori present (authentic standards /
     spiked mix), one per line as "RT,compound_name". Coverage against THIS is
     true accuracy because the truth is independent of both pipeline & manual.

    python eval/regression.py
    python eval/regression.py --truth-set standards.csv

Baselines (vs xt_filled.xlsx, regression-only): top-1>=20, top-5>=23, top-8>=26.
"""
import sys
import argparse
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_DIR))
sys.path.insert(0, str(PROJECT_DIR / "eval"))

REF_MZML = PROJECT_DIR / "output" / "mzml" / "xt6.26.mzML"
CANDIDATES = PROJECT_DIR / "output" / "xt6.26_candidates.json"

# Regression-only floors (agreement with manual top-1; NOT accuracy)
MIN_TOP1 = 20
MIN_TOP5 = 23
MIN_TOP8 = 26

# Dual-evidence verified anchors (strong MS + RI on this column's calibration
# line). Must all stay found in top-8 — a must-not-regress floor.
VERIFIED_ANCHORS = [
    "2-Pentylfuran", "1-Octen-3-one", "(E)-2-Heptenal", "2,3-Octanedione",
    "2,4,6-Trimethylpyridine", "Nonanal", "(E)-2-Octenal", "1-Octen-3-ol",
    "2-Ethyl-1-hexanol", "Benzaldehyde", "(E)-2-Nonenal", "(E)-2-Nonen-1-ol",
    "2-(2-Butoxyethoxy)ethanol", "2,4-Decadienal", "Butylated hydroxytoluene",
]


def _load_truth_set(path):
    """One compound per line: 'RT,name' (RT optional). Lines starting with # skip."""
    rows = []
    for ln in Path(path).read_text(encoding='utf-8').splitlines():
        ln = ln.strip()
        if not ln or ln.startswith('#'):
            continue
        if ',' in ln:
            rt, name = ln.split(',', 1)
            try:
                rows.append({'rt': float(rt), 'name': name.strip(), 'area_pct': 0})
            except ValueError:
                rows.append({'rt': None, 'name': ln, 'area_pct': 0})
        else:
            rows.append({'rt': None, 'name': ln, 'area_pct': 0})
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--truth-set', default=None,
                    help='file of a-priori-known compounds (RT,name) for TRUE accuracy')
    args = ap.parse_args()

    # 1. Run the pipeline (regenerates candidates JSON)
    from pipeline import run_gcms_pipeline
    run_gcms_pipeline(mzml_files=[str(REF_MZML)], config={'use_nist': True})

    from score import (load_ground_truth, load_candidates, score,
                       DEFAULT_GROUND_TRUTH, names_match)
    user = load_ground_truth(DEFAULT_GROUND_TRUTH)
    peaks = load_candidates(str(CANDIDATES))
    n = len(user)
    top1 = sum(1 for r in score(user, peaks, top_n=1) if r['top1'])
    top5 = sum(1 for r in score(user, peaks, top_n=5) if r['top5'])
    res8 = score(user, peaks, top_n=8)
    top8 = sum(1 for r in res8 if r['top5'])

    # 2. Verified anchors must all be found in top-8
    missing = []
    for a in VERIFIED_ANCHORS:
        if not any(names_match(a, r['user']['name']) and r['top5'] for r in res8):
            missing.append(a)

    print(f"\n{'='*56}")
    print(f"  AGREEMENT vs manual top-1  (REGRESSION-ONLY, NOT accuracy)")
    print(f"    top-1: {top1}/{n}  (floor {MIN_TOP1})")
    print(f"    top-5: {top5}/{n}  (floor {MIN_TOP5})")
    print(f"    top-8: {top8}/{n}  (floor {MIN_TOP8})")
    print(f"  VERIFIED ANCHORS (15, must stay in top-8): "
          f"{len(VERIFIED_ANCHORS)-len(missing)}/{len(VERIFIED_ANCHORS)}")
    for m in missing:
        print(f"    LOST: {m}")

    ok = (top1 >= MIN_TOP1 and top5 >= MIN_TOP5 and top8 >= MIN_TOP8
          and not missing)

    # 3. Optional TRUE accuracy vs an independent truth set (standards/spike)
    if args.truth_set:
        truth = _load_truth_set(args.truth_set)
        # if a truth compound has no RT, match by name anywhere in candidates
        def covered(t, depth):
            if t['rt'] is not None:
                rr = score([t], peaks, top_n=depth)
                return rr[0]['top5']
            return any(names_match(t['name'], c['name'])
                       for pk in peaks for c in pk['candidates'][:depth])
        a5 = sum(1 for t in truth if covered(t, 5))
        a8 = sum(1 for t in truth if covered(t, 8))
        print(f"\n  ACCURACY vs truth set '{Path(args.truth_set).name}'  (REAL accuracy)")
        print(f"    top-5: {a5}/{len(truth)} = {a5*100//max(len(truth),1)}%")
        print(f"    top-8: {a8}/{len(truth)} = {a8*100//max(len(truth),1)}%")

    print(f"  {'PASS' if ok else 'FAIL — regression or lost anchor!'}")
    print(f"{'='*56}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
