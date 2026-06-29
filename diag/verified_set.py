#!/usr/bin/env python3
"""
diag/verified_set.py — build a dual-evidence (MS + RI) confirmed set.
=====================================================================
Uses the zgwt-corrected n-alkane RI ladder (carbon numbers verified by M+).
A user compound is "confirmed" when TWO orthogonal evidences agree:
  (1) spectral: the pipeline's NIST search found that compound at this peak
      with a strong reverse match (RMF >= RMF_MIN), AND
  (2) retention: measured RI (corrected ladder) matches a LITERATURE DB-WAX RI
      for that compound within RI_TOL.
This is the recognized confident-identification standard short of running an
authentic standard of each compound.
"""
import sys, json
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "eval"))
from score import names_match, load_ground_truth, DEFAULT_GROUND_TRUTH, DEFAULT_CANDIDATES, load_candidates

# zgwt-corrected ladder: (RT, RI=100*n), carbon numbers verified by M+ = 14n+2
LADDER = [(7.20,1100),(9.86,1200),(12.56,1300),(15.18,1400),(17.69,1500),
          (20.08,1600),(22.36,1700),(24.54,1800),(26.79,1900),(29.03,2000),
          (31.12,2100),(33.08,2200),(34.99,2300),(37.32,2400),(39.43,2500),
          (41.09,2600),(42.50,2700)]
RT_ARR = np.array([r for r, _ in LADDER]); RI_ARR = np.array([i for _, i in LADDER])

RMF_MIN = 800
RI_TOL = 20

CONTAM = ('silox', 'silan', 'phthal', 'siloxane', 'tert-butyl', 'bht')


def meas_ri(rt):
    return float(np.interp(rt, RT_ARR, RI_ARR))


def load_lit_wax():
    """Literature DB-WAX RI only (ri_dual_column 'wax'), not predicted."""
    p = Path(__file__).resolve().parent.parent / "library" / "ri_dual_column.json"
    d = json.load(open(p, encoding='utf-8'))
    return {k.lower(): v for k, v in d.get('wax', {}).items()}


def lit_ri_for(name, wax):
    n = name.lower().strip()
    if n in wax:
        return wax[n]
    for k, v in wax.items():       # fall back to name-equivalence
        if names_match(name, k):
            return v
    return None


def main():
    user = load_ground_truth(DEFAULT_GROUND_TRUTH)
    peaks = load_candidates(str(DEFAULT_CANDIDATES))
    wax = load_lit_wax()

    def peak_at(rt, tol=0.3):
        best, bd = None, tol
        for pk in peaks:
            d = abs(pk['rt'] - rt)
            if d < bd:
                bd, best = d, pk
        return best

    confirmed, ms_only, ri_only, neither = [], [], [], []
    for up in user:
        pk = peak_at(up['rt'])
        mri = meas_ri(up['rt'])
        # spectral evidence: best candidate matching the user compound
        rmf, litri = 0, None
        if pk:
            for c in pk['candidates']:
                if names_match(up['name'], c['name']):
                    rmf = max(rmf, c.get('rmf', 0))
                    litri = litri or lit_ri_for(c['name'], wax)
        if litri is None:
            litri = lit_ri_for(up['name'], wax)
        ms_ok = rmf >= RMF_MIN
        ri_ok = litri is not None and abs(mri - litri) <= RI_TOL
        contam = any(c in up['name'].lower() for c in CONTAM)
        rec = {'name': up['name'], 'rt': up['rt'], 'measRI': round(mri),
               'litRI': litri, 'rmf': rmf,
               'dRI': (round(mri - litri) if litri else None)}
        if ms_ok and ri_ok and not contam:
            confirmed.append(rec)
        elif ms_ok:
            ms_only.append(rec)
        elif ri_ok:
            ri_only.append(rec)
        else:
            neither.append(rec)

    print("=" * 72)
    print(f"  DUAL-EVIDENCE CONFIRMED  (RMF>={RMF_MIN} AND |measRI-litRI|<={RI_TOL})")
    print("=" * 72)
    for r in sorted(confirmed, key=lambda x: x['rt']):
        print(f"  RT={r['rt']:5.2f}  {r['name']:<34} RMF={r['rmf']:<4} "
              f"measRI={r['measRI']} litRI={r['litRI']} dRI={r['dRI']:+d}")
    print(f"\n  CONFIRMED: {len(confirmed)}/{len(user)}")

    print("\n" + "=" * 72)
    print("  MS-STRONG, RI not confirming (split: no-lit-RI vs RI-CONFLICT)")
    print("=" * 72)
    for r in sorted(ms_only, key=lambda x: x['rt']):
        if r['litRI'] is None:
            tag = "no lit RI"
        else:
            tag = f"RI CONFLICT dRI={r['dRI']:+d} (litRI={r['litRI']})"
        print(f"  RT={r['rt']:5.2f}  {r['name']:<34} RMF={r['rmf']:<4} measRI={r['measRI']}  {tag}")

    print(f"\n  CONFIRMED {len(confirmed)} | MS-strong-noRI/conflict {len(ms_only)} | "
          f"RI-ok-weakMS {len(ri_only)} | neither {len(neither)}")
    print("\n  (MS-only = spectral strong, but compound lacks a literature DB-WAX RI")
    print("   or RI disagrees — includes the branched alkanes & E/Z isomers.)")


if __name__ == "__main__":
    main()
