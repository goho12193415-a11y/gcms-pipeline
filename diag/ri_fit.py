#!/usr/bin/env python3
"""
diag/ri_fit.py — is measured RI (zgwt ladder) vs literature RI a clean line?
============================================================================
Fits literature_RI = a*measured_RI + b on compounds whose identity is certain
from their distinctive EI spectra (high RMF). If R^2 is high, the column<->lit
transform is well-defined and RI outliers are meaningful (-> verified set by
residual). If R^2 is low, the literature table is too dirty for RI verification.
Robust: fits, drops gross outliers (bad lit entries), refits.
"""
import sys, json
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "eval"))
from score import names_match, load_ground_truth, DEFAULT_GROUND_TRUTH, DEFAULT_CANDIDATES, load_candidates

LADDER = [(7.20,1100),(9.86,1200),(12.56,1300),(15.18,1400),(17.69,1500),
          (20.08,1600),(22.36,1700),(24.54,1800),(26.79,1900),(29.03,2000),
          (31.12,2100),(33.08,2200),(34.99,2300),(37.32,2400),(39.43,2500),
          (41.09,2600),(42.50,2700)]
RT_ARR = np.array([r for r,_ in LADDER]); RI_ARR = np.array([i for _,i in LADDER])
RMF_MIN = 850


def meas_ri(rt): return float(np.interp(rt, RT_ARR, RI_ARR))


def load_lit_wax():
    p = Path(__file__).resolve().parent.parent / "library" / "ri_dual_column.json"
    d = json.load(open(p, encoding='utf-8'))
    return {k.lower(): v for k, v in d.get('wax', {}).items()}


def lit_ri_for(name, wax):
    n = name.lower().strip()
    if n in wax: return wax[n]
    for k, v in wax.items():
        if names_match(name, k): return v
    return None


def main():
    user = load_ground_truth(DEFAULT_GROUND_TRUTH)
    peaks = load_candidates(str(DEFAULT_CANDIDATES))
    wax = load_lit_wax()

    def peak_at(rt, tol=0.3):
        best, bd = None, tol
        for pk in peaks:
            d = abs(pk['rt']-rt)
            if d < bd: bd, best = d, pk
        return best

    pts = []   # (name, rt, measRI, litRI, rmf)
    for up in user:
        pk = peak_at(up['rt'])
        rmf, litri = 0, None
        if pk:
            for c in pk['candidates']:
                if names_match(up['name'], c['name']):
                    rmf = max(rmf, c.get('rmf',0))
                    litri = litri or lit_ri_for(c['name'], wax)
        litri = litri or lit_ri_for(up['name'], wax)
        if litri and rmf >= RMF_MIN:
            pts.append([up['name'], up['rt'], meas_ri(up['rt']), litri, rmf])

    x = np.array([p[2] for p in pts], float)  # measRI
    y = np.array([p[3] for p in pts], float)  # litRI

    def fit(x, y):
        a, b = np.polyfit(x, y, 1)
        pred = a*x + b
        ss = 1 - np.sum((y-pred)**2)/np.sum((y-np.mean(y))**2)
        return a, b, pred, ss

    a, b, pred, r2 = fit(x, y)
    resid = y - pred
    # robust: drop |resid|>60 (bad lit entries), refit
    keep = np.abs(resid) <= 60
    a2, b2, pred2, r2b = fit(x[keep], y[keep])
    resid_all = y - (a2*x + b2)

    print(f"n points (RMF>={RMF_MIN}, has litRI): {len(pts)}")
    print(f"raw fit:    litRI = {a:.3f}*measRI + {b:.0f}   R2={r2:.3f}")
    print(f"robust fit: litRI = {a2:.3f}*measRI + {b2:.0f}   R2={r2b:.3f}  (dropped {int((~keep).sum())} outliers)")
    print("=" * 78)
    print(f"  {'compound':<34}{'measRI':>7}{'litRI':>7}{'resid':>7}  flag")
    print("=" * 78)
    for (nm, rt, mri, lri, rmf), rr in sorted(zip(pts, resid_all), key=lambda z: z[0][1]):
        flag = "OK" if abs(rr) <= 30 else ("susp" if abs(rr) <= 60 else "OUTLIER")
        print(f"  {nm:<34}{mri:7.0f}{lri:7.0f}{rr:+7.0f}  {flag}")


if __name__ == "__main__":
    main()
