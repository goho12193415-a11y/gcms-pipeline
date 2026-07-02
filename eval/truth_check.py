#!/usr/bin/env python3
"""
eval/truth_check.py — REAL accuracy against authentic standards.
================================================================
Unlike eval/score.py / regression.py (which measure *agreement with the
fallible manual picks* on the xt6.26 reference), this scores the pipeline's
candidates for ANY sample against an independent truth set — the compounds you
KNOW are present because you ran an authentic standard / spiked mix. That makes
the resulting top-1/5/8 numbers a true accuracy, not a consistency ratio.

Two-step workflow:
  1) Run the standard through the pipeline to get its candidates JSON:
       python pipeline.py --files STD.raw --standard alkanes.raw --nist --output output
     -> output/STD_candidates.json
  2) Score it against your truth list:
       python eval/truth_check.py --candidates output/STD_candidates.json --truth truth.csv

Truth file: one compound per line, "RT,name" (RT in minutes, optional).
  - With RT   -> we locate the nearest pipeline peak (within --rt-tol) and check
                 at what rank the compound appears there.
  - Without RT -> we search every peak and report the best (lowest) rank found
                 anywhere. Weaker: a name may coincidentally match an unrelated
                 peak, so give RTs when you know them.
  # lines starting with '#' are ignored.

    Nonanal
    9.85,Benzaldehyde
    12.30,Limonene
"""
import sys
import argparse
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_DIR))
sys.path.insert(0, str(PROJECT_DIR / "eval"))

# Emit UTF-8 so the Chinese labels don't crash on a legacy (GBK) Windows console.
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

from score import names_match, load_candidates  # reuse the validated matcher


def load_truth(path):
    """One compound per line: 'RT,name' (RT optional). '#' comments skipped."""
    rows = []
    for ln in Path(path).read_text(encoding='utf-8').splitlines():
        ln = ln.strip()
        if not ln or ln.startswith('#'):
            continue
        if ',' in ln:
            rt, name = ln.split(',', 1)
            try:
                rows.append({'rt': float(rt), 'name': name.strip()})
                continue
            except ValueError:
                pass  # not "RT,name" — treat the whole line as a name
        rows.append({'rt': None, 'name': ln})
    return rows


def find_rank(truth_name, cands):
    """1-based rank at which truth_name matches in this candidate list, else None."""
    for i, c in enumerate(cands):
        if names_match(truth_name, c.get('name', '')):
            return i + 1
    return None


def check(truth, peaks, rt_tol):
    """For each truth compound resolve (rank, peak) it was found at."""
    out = []
    for t in truth:
        if t['rt'] is not None:
            # nearest pipeline peak by RT
            best, best_d = None, rt_tol
            for pk in peaks:
                d = abs(pk['rt'] - t['rt'])
                if d < best_d:
                    best_d, best = d, pk
            if best is None:
                out.append({'truth': t, 'rank': None, 'peak_rt': None,
                            'drt': None, 'rank1_name': None, 'reason': 'no peak near RT'})
                continue
            cands = best.get('candidates', [])
            rank = find_rank(t['name'], cands)
            r1 = cands[0]['name'] if cands else None
            reason = None
            if rank is None:
                reason = (f"not in top-8 at peak {best['rt']} "
                          f"(rank1 there: {r1})" if r1 else "peak has no candidates")
            out.append({'truth': t, 'rank': rank, 'peak_rt': best['rt'],
                        'drt': round(best_d, 3), 'rank1_name': r1,
                        'reason': reason})
        else:
            # no RT: best rank anywhere
            best_rank, best_pk = None, None
            for pk in peaks:
                r = find_rank(t['name'], pk.get('candidates', []))
                if r is not None and (best_rank is None or r < best_rank):
                    best_rank, best_pk = r, pk
            out.append({'truth': t, 'rank': best_rank,
                        'peak_rt': best_pk['rt'] if best_pk else None,
                        'drt': None, 'rank1_name': None,
                        'reason': None if best_rank else 'not found in any peak'})
    return out


def main():
    ap = argparse.ArgumentParser(description='REAL accuracy vs authentic standards')
    ap.add_argument('--candidates', required=True,
                    help='pipeline output/<sample>_candidates.json')
    ap.add_argument('--truth', required=True,
                    help='truth list, one compound per line "RT,name" (RT optional)')
    ap.add_argument('--rt-tol', type=float, default=0.30,
                    help='RT match tolerance in minutes (default 0.30)')
    args = ap.parse_args()

    peaks = load_candidates(args.candidates)
    truth = load_truth(args.truth)
    res = check(truth, peaks, args.rt_tol)
    n = len(truth)

    def acc(depth):
        return sum(1 for r in res if r['rank'] is not None and r['rank'] <= depth)

    t1, t5, t8 = acc(1), acc(5), acc(8)
    pct = lambda x: f"{x*100//max(n,1)}%"

    print(f"\n{'='*60}")
    print(f"  TRUTH-SET ACCURACY  (真实准确率 vs 已知标准品)")
    print(f"  candidates: {Path(args.candidates).name}")
    print(f"  truth:      {Path(args.truth).name}   ({n} compounds)")
    print(f"  {'-'*56}")
    print(f"    top-1: {t1}/{n} = {pct(t1)}   (管道首选就对)")
    print(f"    top-5: {t5}/{n} = {pct(t5)}")
    print(f"    top-8: {t8}/{n} = {pct(t8)}")
    print(f"  {'-'*56}")
    print(f"  逐条明细:")
    for r in res:
        t = r['truth']
        tag = f"RT{t['rt']}" if t['rt'] is not None else "no-RT"
        if r['rank'] is None:
            print(f"    MISS      {t['name'][:34]:34s} [{tag}] - {r['reason']}")
        else:
            loc = (f"peak {r['peak_rt']} drt{r['drt']}" if t['rt'] is not None
                   else f"@peak {r['peak_rt']}")
            extra = ''
            if t['rt'] is not None and r['rank'] > 1 and r['rank1_name']:
                extra = f"  (rank1 was: {r['rank1_name'][:28]})"
            print(f"    hit@{r['rank']:<4}  {t['name'][:34]:34s} [{tag}] {loc}{extra}")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()
