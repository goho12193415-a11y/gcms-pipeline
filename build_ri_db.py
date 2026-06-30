#!/usr/bin/env python3
"""
build_ri_db.py — (re)build the local RI database from the user's NIST ri.dat.
============================================================================
Parses the NIST RI library's ri.dat (the user's own licensed data; output JSON
is local-only and git-ignored, NOT redistributed). The previous database kept
only phases labelled exactly "DB-Wax"; this also folds in NIST's broad
"Standard polar" column class plus all equivalent polar phases (Carbowax 20M,
Supelcowax, HP-INNOWax, FFAP, PEG ...), and likewise for non-polar (DB-5),
substantially increasing coverage. Per compound it stores the MEDIAN RI.

    python build_ri_db.py "<...>/nist_ri/ri.dat" [out.json]

ri.dat record layout (text fields in sequence): C-ID, name, formula, then per
measurement: RI, "Capillary"/"Packed", column-class, specific-phase, length...
"""
import sys, re, json, statistics
from pathlib import Path

POLAR_CLASS = {'standard polar'}
NONPOLAR_CLASS = {'standard non-polar', 'semi-standard non-polar'}
POLAR_PH = ('carbowax', 'db-wax', 'supelcowax', 'innowax', 'stabilwax', 'ffap',
            'polyethylene glycol', 'cp-wax', 'at-wax', 'zb-wax', 'hp-20m', 'peg',
            'omegawax', 'solgel-wax', 'nukol')
NONPOLAR_PH = ('db-5', 'hp-5', 'se-30', 'ov-101', 'db-1', 'hp-1', 'ov-1', 'se-54',
               'cp-sil 5', 'cp-sil 8', 'rtx-5', 'zb-5', 'db-5ms', 'hp-5ms', 'se-52',
               'rtx-1', 'mdn-5', 'ultra-1', 'ultra-2', 'pte-5', 'apiezon')


def classify(cls, phase):
    """-> 'wax' / 'db5' / None, using NIST column-class then specific phase."""
    c, p = cls.lower().strip(), phase.lower().strip()
    if c in POLAR_CLASS:
        return 'wax'
    if c in NONPOLAR_CLASS:
        return 'db5'
    if any(x in p for x in POLAR_PH):
        return 'wax'
    if any(x in p for x in NONPOLAR_PH):
        return 'db5'
    return None


def is_cid(s):
    return len(s) > 2 and s[0] == 'C' and s[1:].isdigit()


def is_ri(s):
    try:
        f = float(s); return 350.0 <= f <= 4600.0
    except ValueError:
        return False


def main():
    ri_dat = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) > 2 else 'library/ri_dual_column.json'
    data = open(ri_dat, 'rb').read()
    S = [m.decode('latin1') for m in re.findall(rb'[ -~]{2,}', data)]
    print(f"  {len(S):,} strings extracted from {Path(ri_dat).name}")

    wax, db5 = {}, {}          # name_lower -> list of RI
    cur = None
    N = len(S)
    for i, s in enumerate(S):
        if is_cid(s):
            cur = S[i + 1].strip() if i + 1 < N else None
            continue
        if cur and is_ri(s):
            ri = float(s)
            # measurement fields follow the RI; find class + specific phase nearby
            window = S[i + 1:i + 5]
            cls = next((w for w in window
                        if w.lower().strip() in POLAR_CLASS | NONPOLAR_CLASS), '')
            phase = ''
            for w in window:
                if any(x in w.lower() for x in POLAR_PH + NONPOLAR_PH):
                    phase = w; break
            kind = classify(cls, phase)
            if kind == 'wax':
                wax.setdefault(cur.lower(), []).append(ri)
            elif kind == 'db5':
                db5.setdefault(cur.lower(), []).append(ri)

    wax_med = {k: int(round(statistics.median(v))) for k, v in wax.items()}
    db5_med = {k: int(round(statistics.median(v))) for k, v in db5.items()}
    print(f"  WAX (polar):     {len(wax_med):,} compounds")
    print(f"  DB-5 (non-polar):{len(db5_med):,} compounds")
    for nm in ('nonanal', 'octanal', 'benzaldehyde', '1-octen-3-ol', 'hexanal'):
        print(f"    {nm:<16} wax={wax_med.get(nm)} db5={db5_med.get(nm)}")
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    json.dump({'wax': wax_med, 'db5': db5_med}, open(out, 'w', encoding='utf-8'),
              ensure_ascii=False)
    print(f"  wrote {out}")


if __name__ == "__main__":
    main()
