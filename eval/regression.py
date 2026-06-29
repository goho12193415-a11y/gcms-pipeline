#!/usr/bin/env python3
"""
eval/regression.py — coverage regression guard.
===============================================
Runs the pipeline on the reference sample (xt6.26) and asserts the top-5
coverage has not dropped below the locked baseline. Exits non-zero on
regression so it can gate any future change.

    python eval/regression.py

Baselines (eval/score.py vs xt_filled.xlsx, 49 manual IDs):
    2026-06-29 air-strip + name fixes:        top-1 >= 20, top-5 >= 23
    2026-06-29 + replib merge (top-8 review): top-8 >= 26
"""
import sys
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_DIR))
sys.path.insert(0, str(PROJECT_DIR / "eval"))

REF_MZML = PROJECT_DIR / "output" / "mzml" / "xt6.26.mzML"
CANDIDATES = PROJECT_DIR / "output" / "xt6.26_candidates.json"

MIN_TOP1 = 20
MIN_TOP5 = 23
MIN_TOP8 = 26


def main():
    # 1. Run the pipeline (regenerates the candidates JSON)
    from pipeline import run_gcms_pipeline
    run_gcms_pipeline(mzml_files=[str(REF_MZML)], config={'use_nist': True})

    # 2. Score
    from score import load_ground_truth, load_candidates, score, DEFAULT_GROUND_TRUTH
    user = load_ground_truth(DEFAULT_GROUND_TRUTH)
    peaks = load_candidates(str(CANDIDATES))
    n = len(user)
    top1 = sum(1 for r in score(user, peaks, top_n=1) if r['top1'])
    top5 = sum(1 for r in score(user, peaks, top_n=5) if r['top5'])
    top8 = sum(1 for r in score(user, peaks, top_n=8) if r['top5'])

    print(f"\n{'='*50}")
    print(f"  REGRESSION CHECK")
    print(f"  top-1: {top1}/{n}  (floor {MIN_TOP1})")
    print(f"  top-5: {top5}/{n}  (floor {MIN_TOP5})")
    print(f"  top-8: {top8}/{n}  (floor {MIN_TOP8})")
    ok = top1 >= MIN_TOP1 and top5 >= MIN_TOP5 and top8 >= MIN_TOP8
    print(f"  {'PASS' if ok else 'FAIL — coverage regressed!'}")
    print(f"{'='*50}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
