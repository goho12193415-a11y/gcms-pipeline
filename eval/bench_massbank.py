#!/usr/bin/env python3
"""
eval/bench_massbank.py — independent identification accuracy vs MassBank EI.
===========================================================================
An accuracy check that needs NO lab work and NO analyst picks: it feeds
INDEPENDENTLY-acquired EI mass spectra of KNOWN compounds (MassBank of North
America / MassBank EU) through this project's NIST search, and scores how often
the correct compound is returned — with identity judged by STRUCTURE (CAS from
PubChem, keyed on the record's InChIKey), not by name string. That removes the
naming-artifact undercount that plagues name-only scoring
(e.g. "m-hydroxy nitrobenzene" == NIST "Phenol, 3-nitro-").

Why this is meaningful (and its limits):
  + independent of NIST (MassBank spectra are not NIST's own) -> not circular
  + structure-based identity (CAS/InChIKey) -> no name-matching undercount
  - single-compound LIBRARY spectra are cleaner than real chromatographic peaks
    (co-elution/noise), so this is an UPPER bound vs real samples
  - MassBank's EI set skews to derivatized/exotic metabolomics standards; use
    --max-mw / the underivatized filter to approximate real volatile profiling
  - MS-only: it does NOT use retention index. The real pipeline's RI cross-check
    disambiguates isomers this test cannot, so real dual-evidence accuracy on
    RI-supported compounds is HIGHER than the number here.

Headline result (2026-07, NIST14 mainlib, underivatized, MW<=200, n=300):
  raw:                       top-1 ~60%, top-5 ~70%
  given compound retrievable: retrievable 71%, top-1 89%, top-5 ~100%
  => the search RANKING is near-ceiling; the gap is library coverage + EI
     isomer ambiguity, not the matching algorithm.

Data is downloaded on demand and cached locally; none of it is committed
(MassBank/NIST-derived — see .gitignore). Requires network for the first run
(MassBank release asset + PubChem REST) and the NIST DLL (via nist_engine).

Usage:
  python eval/bench_massbank.py                 # defaults: underivatized, MW<=200, n=300
  python eval/bench_massbank.py --n 1000 --max-mw 250
  python eval/bench_massbank.py --msp path/to/MassBank_NISTformat.msp
"""
import os, sys, re, json, time, random, argparse, urllib.request
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_DIR))
sys.path.insert(0, str(PROJECT_DIR / "eval"))

CACHE_DIR = PROJECT_DIR / "eval" / "_bench_cache"   # gitignored
MSP_URL = ("https://github.com/MassBank/MassBank-data/releases/download/"
           "2026.03/MassBank_NISTformat.msp")
DERIV_MARKERS = ('tms', 'silyl', 'silox', 'silane', 'trifluoroacet', 'pentafluoro',
                 'heptafluoro', 'tbs', 'tert-butyldimethyl', 'derivative', ' deriv',
                 'fluoroacetyl', 'pfp', 'pfb', 'trimethylsilyl',
                 '4tms', '3tms', '2tms', '5tms')
_AW = {'C': 12, 'H': 1, 'O': 16, 'N': 14, 'S': 32, 'Cl': 35.5, 'Br': 80,
       'F': 19, 'P': 31, 'Si': 28, 'I': 127}


def _mw(formula):
    t = 0
    for el, ct in re.findall(r'([A-Z][a-z]?)(\d*)', formula or ''):
        if el in _AW:
            t += _AW[el] * (int(ct) if ct else 1)
    return t or None


def _ncas(c):
    """Normalize a CAS number (strip leading zeros per segment) or None."""
    if not c or c == '0-0-0':
        return None
    p = c.strip().split('-')
    if len(p) != 3 or not all(x.isdigit() for x in p):
        return None
    return '-'.join(str(int(x)) for x in p)


def _canon(name):
    """Canonical form that inverts NIST's 'parent, substituent-' naming so it
    lines up with assembled names, then strips to alphanumerics."""
    s = (name or '').lower().strip()
    s = re.sub(r'\.(alpha|beta|gamma|delta|omega)\.', '', s)
    for g in 'αβγδ':
        s = s.replace(g, '')
    if ', ' in s:
        parent, rest = s.split(', ', 1)
        s = rest + parent
    s = re.sub(r'\(.*?\)', '', s)
    return re.sub(r'[^a-z0-9]', '', s)


def download_msp(dest):
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"[bench] downloading MassBank release ({MSP_URL.split('/')[-1]}, ~134MB)...")
    urllib.request.urlretrieve(MSP_URL, dest)
    return dest


def parse_ei(msp_path):
    """Parse an NIST-format MSP, keeping only EI (GC-EI) records with a name."""
    text = Path(msp_path).read_text(encoding='utf-8', errors='ignore')
    kept, cur, pk, reading = [], {}, [], False
    def close():
        it = cur.get('instrument_type', '')
        if 'EI' in it and 'ESI' not in it and cur.get('name') and pk:
            kept.append({'name': cur['name'], 'inchikey': cur.get('inchikey', ''),
                         'formula': cur.get('formula', ''), 'peaks': pk[:]})
    for line in text.split('\n'):
        line = line.rstrip('\r')
        if line.strip() == '':
            close(); cur, pk, reading = {}, [], False; continue
        if reading:
            for tok in line.split(';'):
                p = tok.split()
                if len(p) >= 2:
                    try: pk.append((float(p[0]), float(p[1])))
                    except ValueError: pass
        elif ':' in line:
            k, v = line.split(':', 1)
            k = k.strip().lower().replace(' ', '_'); cur[k] = v.strip()
            if k == 'num_peaks': reading = True
    close()
    return kept


class PubChemCAS:
    """InChIKey -> (set of CAS, set of synonym names), cached to disk."""
    def __init__(self, path):
        self.path = Path(path)
        self.cache = json.loads(self.path.read_text()) if self.path.exists() else {}
    def save(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.path.write_text(json.dumps(self.cache))
    def get(self, ik):
        if ik in self.cache:
            cas, nm = self.cache[ik]; return set(cas), set(nm)
        u = (f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/"
             f"{ik}/synonyms/JSON")
        try:
            r = json.load(urllib.request.urlopen(
                urllib.request.Request(u, headers={'User-Agent': 'gcms-bench'}), timeout=20))
            syn = r['InformationList']['Information'][0].get('Synonym', [])
        except Exception:
            syn = []
        cas, nm = set(), set()
        for x in syn:
            n = _ncas(x)
            (cas.add(n) if n else nm.add(x.lower()))
        self.cache[ik] = [sorted(cas), sorted(nm)]
        time.sleep(0.12)
        return cas, nm


def main():
    ap = argparse.ArgumentParser(description="Independent ID accuracy vs MassBank EI")
    ap.add_argument('--msp', default=str(CACHE_DIR / "MassBank_NISTformat.msp"))
    ap.add_argument('--n', type=int, default=300, help="sample size")
    ap.add_argument('--max-mw', type=float, default=200,
                    help="keep compounds with MW <= this (0 = no cap)")
    ap.add_argument('--keep-derivatized', action='store_true',
                    help="do NOT filter out TMS/TFA/etc. derivatives")
    ap.add_argument('--retr-depth', type=int, default=30,
                    help="a compound counts as 'retrievable' if correct is in top-N")
    ap.add_argument('--seed', type=int, default=0)
    args = ap.parse_args()

    from score import names_match
    from nist_engine import NISTSearchEngine
    from pyms.Spectrum import MassSpectrum

    # --- silence the NIST DLL's fd-level stdout spam; keep a clean out stream ---
    out = os.fdopen(os.dup(1), 'w', encoding='utf-8')
    devnull = os.open(os.devnull, os.O_WRONLY); os.dup2(devnull, 1)
    def log(*a): print(*a, file=out); out.flush()

    msp = Path(args.msp)
    if not msp.exists():
        download_msp(msp)
    ei_cache = CACHE_DIR / "massbank_ei.json"
    if ei_cache.exists():
        data = json.loads(ei_cache.read_text())
    else:
        data = parse_ei(msp)
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        ei_cache.write_text(json.dumps(data))
    log(f"[bench] EI records parsed: {len(data)}")

    def keep(d):
        if not d.get('inchikey'):
            return False
        if not args.keep_derivatized and any(b in d['name'].lower() for b in DERIV_MARKERS):
            return False
        if args.max_mw and (_mw(d.get('formula', '')) or 9999) > args.max_mw:
            return False
        return True
    pool = [d for d in data if keep(d)]
    random.seed(args.seed); random.shuffle(pool)
    sample = pool[:args.n]
    log(f"[bench] filtered pool={len(pool)}  scoring n={len(sample)} "
        f"(underivatized={not args.keep_derivatized}, max_mw={args.max_mw})")

    eng = NISTSearchEngine()
    pc = PubChemCAS(CACHE_DIR / "pubchem_cas.json")

    def search(peaks, tn):
        d = {}
        for mz, it in peaks:
            m = int(round(mz))
            if m >= 1: d[m] = d.get(m, 0) + it
        if not d: return []
        mx = max(d.values()); mzs = sorted(d)
        ms = MassSpectrum([float(m) for m in mzs], [int(d[m] / mx * 999) for m in mzs])
        try:
            return list(eng.engine.full_search_with_ref_data(ms))[:tn]
        except Exception:
            return []

    def correct_rank(s, hits):
        pcas, pnames = pc.get(s['inchikey'])
        mb = pnames | {s['name'].lower()}; cmb = _canon(s['name'])
        for i, h in enumerate(hits):
            ref = h[1]; hc = _ncas(getattr(ref, 'cas', ''))
            names = [getattr(ref, 'name', '')] + list(getattr(ref, 'synonyms', None) or [])
            if hc and hc in pcas:
                return i + 1
            for x in names:
                if x.lower() in mb or names_match(s['name'], x) or (cmb and cmb == _canon(x)):
                    return i + 1
        return None

    n = t1 = t5 = t8 = retr = c1 = c5 = c8 = 0
    for s in sample:
        hits = search(s['peaks'], args.retr_depth)
        if not hits:
            continue
        n += 1
        r = correct_rank(s, hits)
        if r is not None:
            if r == 1: t1 += 1
            if r <= 5: t5 += 1
            if r <= 8: t8 += 1
            retr += 1
            if r == 1: c1 += 1
            if r <= 5: c5 += 1
            if r <= 8: c8 += 1
    pc.save()

    pct = lambda a, b: f"{a * 100 // max(b, 1)}%"
    log(f"\n{'='*60}")
    log(f"  MassBank-EI identification accuracy (structure-verified, MS-only)")
    log(f"  scored n={n}   library=NIST via nist_engine")
    log(f"  {'-'*56}")
    log(f"  RAW (all scored):")
    log(f"    top-1 {t1}/{n} = {pct(t1,n)}    top-5 {t5}/{n} = {pct(t5,n)}    "
        f"top-8 {t8}/{n} = {pct(t8,n)}")
    log(f"  RETRIEVABLE (correct compound within top-{args.retr_depth}): "
        f"{retr}/{n} = {pct(retr,n)}")
    log(f"  GIVEN retrievable (search-ranking quality):")
    log(f"    top-1 {c1}/{retr} = {pct(c1,retr)}    top-5 {c5}/{retr} = {pct(c5,retr)}"
        f"    top-8 {c8}/{retr} = {pct(c8,retr)}")
    log(f"  NOTE: library-quality spectra (upper bound vs real peaks); MS-only")
    log(f"        (real pipeline adds RI cross-check -> higher on isomers).")
    log(f"{'='*60}")


if __name__ == "__main__":
    main()
