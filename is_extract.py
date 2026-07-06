#!/usr/bin/env python3
"""
is_extract.py — targeted internal-standard (IS) extraction.
===========================================================
An internal standard that co-elutes with a large peak is missed by ordinary
peak detection. This pulls a chosen ion's chromatogram (EIC), locates and
integrates the IS peak inside a retention-time window, and tabulates the
response per sample — the consistent number you normalise other compounds to.
Useful when a batch has only an internal standard and no calibration standard.

Works on any format the pipeline reads (Thermo .RAW / .mzML / Shimadzu .qgd /
Agilent .D) via step1_parse.load_sample. Default target: 2-octanol (m/z 45).

    python is_extract.py --dir "D:\\...\\folder" --mz 45 --rt-min 11.8 --rt-max 13.2
    python is_extract.py --files a.qgd b.D --mz 45 --out is.xlsx

A confirmation-ion ratio (default m/z 55/45) separates a real IS peak from
noise: genuine 2-octanol sits near 0.4; a wildly different ratio means the tiny
signal in the window is not the standard.
"""
import os, sys, glob, argparse
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def eic(scan_list, target, tol=0.5):
    e = np.zeros(len(scan_list))
    for i, s in enumerate(scan_list):
        m = s['mz']
        if len(m):
            k = (m >= target - tol) & (m <= target + tol)
            if k.any():
                e[i] = s['intensity'][k].sum()
    return e


def integrate_peak(rt, e, lo, hi):
    win = np.where((rt >= lo) & (rt <= hi))[0]
    if len(win) == 0:
        return None
    sub = e[win]
    base = float(np.median(sub))
    apex = int(win[int(np.argmax(sub))])
    thr = base + 0.05 * (e[apex] - base)
    l = apex
    while l - 1 >= win[0] and e[l - 1] > thr:
        l -= 1
    r = apex
    while r + 1 <= win[-1] and e[r + 1] > thr:
        r += 1
    seg = np.clip(e[l:r + 1] - base, 0, None)
    dt = np.diff(rt[l:r + 1])
    area = float(np.sum((seg[:-1] + seg[1:]) / 2 * dt)) if len(seg) > 1 else 0.0
    return {'rt': float(rt[apex]), 'height': float(e[apex] - base),
            'area': area, 'rt_start': float(rt[l]), 'rt_end': float(rt[r])}


def extract_internal_standard(sample_paths, mz=45.0, rt_min=11.8, rt_max=13.2,
                              confirm=55.0, name='2-octanol', log=print):
    """Return a DataFrame of the IS response per sample. `log` receives progress."""
    from step1_parse import load_sample
    rows = []
    for p in sample_paths:
        base = os.path.basename(str(p).rstrip('/\\'))
        nm = base[:-2] if base.lower().endswith('.d') else os.path.splitext(base)[0]
        try:
            d = load_sample(p)
        except Exception as ex:
            log(f"  [IS-FAIL] {nm}: {str(ex)[:120]}")
            continue
        rt = np.asarray(d['rt'], float)
        q = eic(d['scan_list'], mz)
        c = eic(d['scan_list'], confirm)
        pk = integrate_peak(rt, q, rt_min, rt_max)
        if pk is None:
            rows.append({'sample': nm, 'IS_RT': None, 'IS_area': 0, 'IS_height': 0,
                         'confirm_ratio': None, 'is_blank': 'blank' in nm.lower()})
            continue
        ai = int(np.argmin(np.abs(rt - pk['rt'])))
        ratio = round(float(c[ai] / q[ai]), 2) if q[ai] > 0 else None
        rows.append({'sample': nm, 'IS_RT': round(pk['rt'], 3),
                     'IS_area': round(pk['area']), 'IS_height': round(pk['height']),
                     'confirm_ratio': ratio, 'is_blank': 'blank' in nm.lower()})
        log(f"  {nm:26s} RT={pk['rt']:.2f} area={pk['area']:>12.0f} "
            f"m/z{confirm:.0f}/{mz:.0f}={ratio}")
    df = pd.DataFrame(rows)
    smp = df[~df['is_blank']] if 'is_blank' in df else df
    if len(smp) and smp['IS_area'].gt(0).any():
        a = smp['IS_area'][smp['IS_area'] > 0]
        cv = 100 * a.std() / a.mean() if a.mean() else 0
        log(f"[IS] {name}: n={len(a)} 有响应, 均值={a.mean():.0f}, CV={cv:.1f}% "
            f"(CV<~20% 说明内标/进样稳定)")
    return df


def _find_d(folder):
    out = []
    for x in sorted(os.listdir(folder)):
        p = os.path.join(folder, x)
        if x.lower().endswith('.d') and os.path.isfile(os.path.join(p, 'AcqData', 'MSScan.bin')):
            out.append(p)
    return out


def main():
    ap = argparse.ArgumentParser(description="Targeted internal-standard extraction")
    ap.add_argument('--dir', help="folder of .D samples")
    ap.add_argument('--files', nargs='+', help="explicit sample paths")
    ap.add_argument('--mz', type=float, default=45.0)
    ap.add_argument('--rt-min', type=float, default=11.8)
    ap.add_argument('--rt-max', type=float, default=13.2)
    ap.add_argument('--confirm', type=float, default=55.0)
    ap.add_argument('--name', default='2-octanol')
    ap.add_argument('--out', default=None)
    args = ap.parse_args()

    samples = list(args.files or [])
    if args.dir:
        samples += _find_d(args.dir)
    if not samples:
        ap.error("provide --files or --dir")
    print(f"[IS] {len(samples)} samples; m/z {args.mz:.0f} in "
          f"{args.rt_min}-{args.rt_max} min ({args.name})")
    df = extract_internal_standard(samples, args.mz, args.rt_min, args.rt_max,
                                   args.confirm, args.name)
    out = args.out or os.path.join(args.dir or '.', f"internal_standard_{args.name}.xlsx")
    df.to_excel(out, index=False)
    print(f"[IS] wrote {out}")


if __name__ == "__main__":
    main()
