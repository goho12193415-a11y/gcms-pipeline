"""
NIST MS Search Engine v3 — Hybrid RI + MS Search
==================================================
For each peak:
  1. NIST MS search → top-10 candidates by spectral match
  2. RI database search → compounds within +/-50 RI of measured value
  3. Merge both lists — RI-matched compounds added even if MS rank is low
  4. Re-rank: compounds in BOTH lists get top priority
"""
import json, re, numpy as np
from pathlib import Path
from pyms.Spectrum import MassSpectrum
import pyms_nist_search

# Air / background ions that dominate low-abundance peak spectra and wreck the
# NIST match (root cause of the Bm-class failures). Defined in config so it is
# tunable in one place; frozenset here for fast membership tests.
try:
    from config import AIR_IONS as _AIR
    AIR_IONS = frozenset(_AIR)
except ImportError:
    AIR_IONS = frozenset({17, 18, 28, 32, 40, 44})

try:
    from config import COELUTION_RATIO as _COELUTION_RATIO
except ImportError:
    _COELUTION_RATIO = 1.3


def _canon_alkane(name):
    """Canonicalize IUPAC<->common branched-alkane naming to one key.
    'Dodecane, 2,6,10-trimethyl-' / '2,6,10-Trimethyldodecane' -> 'dodecane-2,6,10-trimethyl'"""
    n = name.lower().strip()
    m = re.match(r'([a-z]+),?\s*([\d,]+-[a-z]+)', n)
    if m:
        return f'{m.group(1).strip().rstrip(",")}-{m.group(2).strip().rstrip("-")}'
    m = re.match(r'([\d,]+-[a-z]+)([a-z]+)', n)
    if m:
        return f'{m.group(2)}-{m.group(1)}'
    return n


_ALK_BAD = ('ene', 'yne', 'ol', 'one', 'al,', 'oic', 'acid', 'ester', 'amine',
            'oxy', 'sil', 'phen', 'furan', 'ketone', 'ate', 'bromo', 'chloro',
            'fluoro', 'iodo', 'adamant', 'nitro', 'cyclo', 'thio', 'naphthalen')


def _is_branched_alkane(name):
    """True for methyl-substituted saturated alkanes (no other functional group)."""
    n = name.lower()
    if 'methyl' not in n or 'ane' not in n:
        return False
    return not any(b in n for b in _ALK_BAD)


def _is_alkane_like(name):
    """True for any saturated acyclic alkane (straight or methyl-branched).
    Used to detect that a peak's spectrum reads as an alkane, which is the
    cue that the exact (branched) homolog may be missing from the hits."""
    n = name.lower().strip()
    if not re.search(r'[a-z]ane\b', n):
        return False
    return not any(b in n for b in _ALK_BAD)


class NISTSearchEngine:
    def __init__(self, lib_path=None, work_dir=None, lib_type=None):
        if lib_type is None:
            lib_type = pyms_nist_search.NISTMS_MAIN_LIB
        if lib_path is None:
            candidates = [
                Path(r"C:\Users\go ho\Desktop\MSSEARCH\mainlib"),
                Path(r"C:\Users\go ho\Desktop\GCMS_Software\谱库\mainlib"),
            ]
            lib_path = None
            for c in candidates:
                if c.exists():
                    lib_path = str(c)
                    break
            if lib_path is None:
                lib_path = r"C:\Users\go ho\Desktop\MSSEARCH\mainlib"
        if work_dir is None:
            import tempfile, os
            work_dir = os.path.join(tempfile.gettempdir(), "nist_work")
            os.makedirs(work_dir, exist_ok=True)

        # NOTE: pyms_nist_search wraps a stateful DLL with a single global
        # "active library" — two Engine instances in one process corrupt each
        # other (the second hijacks the DLL). So replib is NOT loaded here;
        # it is searched in a separate process (replib_pass.py) and merged
        # offline by the pipeline. This engine is mainlib-only.
        self.engine = pyms_nist_search.Engine(lib_path, lib_type, work_dir)
        self.lib_type = lib_type
        self._load_ri_database()

    def _load_ri_database(self):
        """Load dual-column RI database."""
        self.ri_db5 = {}   # {name_lower: ri}
        self.ri_wax = {}   # {name_lower: ri}

        from config import PROJECT_DIR as project_dir
        paths = [
            project_dir / "library" / "ri_dual_column.json",
            project_dir / "library" / "ri_nist_full.json",
            project_dir / "library" / "ri_enriched.json",
        ]
        for path in paths:
            try:
                with open(path, 'r', encoding='utf-8') as f:
                    data = json.load(f)
                self.ri_db5.update(data.get('db5', {}))
                self.ri_wax.update(data.get('wax', data.get('name_to_ri', {})))
                # Also handle old CAS-based format
                by_cas = data.get('by_cas', {})
                for cas, entry in by_cas.items():
                    name = entry.get('name', '').lower()
                    if name and 'DB-5' in entry:
                        self.ri_db5[name] = entry['DB-5']
                    if name and 'DB-WAX' in entry:
                        self.ri_wax[name] = entry['DB-WAX']
            except FileNotFoundError:
                continue

        print(f"  [RI] Dual-column: {len(self.ri_db5):,} DB-5 + "
              f"{len(self.ri_wax):,} DB-WAX")
        self._build_alkane_ri_table()

    def _build_alkane_ri_table(self):
        """Predicted-RI table of branched alkanes (for RI-primary injection),
        loaded from ri_nist_full.json name_to_ri (CNN/Goodner predictions)."""
        self.alkane_ri = []        # sorted list of (ri, canon, display_name)
        from config import PROJECT_DIR as project_dir
        try:
            with open(project_dir / "library" / "ri_nist_full.json",
                      'r', encoding='utf-8') as f:
                name_to_ri = json.load(f).get('name_to_ri', {})
        except FileNotFoundError:
            name_to_ri = {}
        seen = {}
        for nm, ri in name_to_ri.items():
            if _is_branched_alkane(nm):
                c = _canon_alkane(nm)
                if c not in seen:        # keep first (predictions are per-name)
                    seen[c] = (ri, c, nm)
        self.alkane_ri = sorted(seen.values())
        print(f"  [RI] Branched-alkane predicted-RI table: {len(self.alkane_ri)}")

    def _inject_ri_alkanes(self, results, measured_ri,
                           tol=25, max_inject=2):
        """If the spectrum confidently reads 'branched alkane' but the exact
        homolog is absent from the spectral hits, inject the predicted-RI
        candidate(s) closest to the measured RI. Inserted after spectral #1
        so the strong spectral hit stays visible; total length preserved."""
        if not results or measured_ri is None or not self.alkane_ri:
            return results
        top3 = results[:3]
        n_alk = sum(1 for r in top3 if _is_alkane_like(r.get('name', '')))
        if n_alk < 2:                    # spectrum not alkane-dominated -> skip
            return results
        have = {_canon_alkane(r.get('name', '')) for r in results}
        cands = sorted(((abs(ri - measured_ri), ri, c, nm)
                        for ri, c, nm in self.alkane_ri
                        if abs(ri - measured_ri) <= tol and c not in have),
                       key=lambda x: x[0])
        if not cands:
            return results
        keep = max(1, len(results) - min(max_inject, len(cands)))
        injected = []
        for d, ri, c, nm in cands[:max_inject]:
            injected.append({
                'name': nm, 'cas': '', 'fmf': 0, 'rmf': 0,
                'ri_wax': ri, 'ri_db5': None,
                'ri_diff': round(d, 1),
                'source': 'RI-predicted',
                'combined_score': float(results[0].get('rmf', 0)),
            })
        merged = results[:1] + injected + results[1:keep]
        for i, r in enumerate(merged):
            r['rank'] = i + 1
        return merged

    def search_raw_spectrum(self, mz_array, int_array, top_n=5,
                              measured_ri=None):
        """NIST search + RI annotation + homolog injection (no re-ranking)."""
        mz_array = np.asarray(mz_array, dtype=float)
        int_array = np.asarray(int_array, dtype=float)
        if len(mz_array) == 0:
            return []

        # Strip air/background ions before normalization, so the real base
        # peak (not water) drives the match.
        keep = np.array([round(m) not in AIR_IONS for m in mz_array])
        if keep.any():
            mz_array = mz_array[keep]
            int_array = int_array[keep]
        if len(mz_array) == 0:
            return []

        max_i = int_array.max()
        if max_i <= 0:
            return []
        int_norm = (int_array / max_i * 999).astype(int)

        # Sort by m/z (NIST DLL requires sorted m/z)
        sort_idx = np.argsort(mz_array)
        mz_sorted = mz_array[sort_idx]
        int_sorted = int_norm[sort_idx]
        if len(mz_sorted) > 200:
            top_idx = np.argsort(int_sorted)[-200:]
            mz_list = mz_sorted[top_idx].tolist()
            int_list = int_sorted[top_idx].tolist()
        else:
            mz_list = mz_sorted.tolist()
            int_list = int_sorted.tolist()

        try:
            ms = MassSpectrum(mz_list, int_list)
            hits = list(self.engine.full_search_with_ref_data(ms))
        except:
            return []

        # ---- RI re-ranking (optional, disabled by default) ----
        # Enable with: USE_RI_RERANK = True
        # Currently off: CAS table coverage too low to improve overall results
        #
        # if measured_ri is not None and sr_list:
        #     try:
        #         from ri_reranker import rerank_pyms_results
        #         reranked = rerank_pyms_results(...)
        #     ...

        # Standard: top-5 from NIST with RI annotation only
        results = []
        for i, hit in enumerate(hits[:top_n]):
            if not isinstance(hit, tuple) or len(hit) < 2:
                continue
            sr, ref = hit[0], hit[1]
            name = getattr(ref, 'name', '') or ''
            cas = getattr(ref, 'cas', '') or ''
            mf = int(getattr(sr, 'match_factor', 0) or 0)
            rmf = int(getattr(sr, 'reverse_match_factor', 0) or 0)

            nm = name.lower()
            ri_wax = self.ri_wax.get(nm)
            ri_db5 = self.ri_db5.get(nm)
            ri_delta = abs(measured_ri - ri_wax) if (ri_wax and measured_ri) else None

            results.append({
                'name': name, 'cas': cas,
                'fmf': mf, 'rmf': rmf,
                'ri_wax': ri_wax, 'ri_db5': ri_db5,
                'ri_diff': ri_delta,
                # RMF (ignores extra ions) >> FMF (penalises them) means the
                # spectrum carries ions the library compound can't explain -> a
                # mixed/co-eluting peak. Threshold from config (NIST scale).
                'coelution_flag': mf > 0 and (rmf / mf) > _COELUTION_RATIO,
                'source': 'MS',
                'combined_score': float(rmf),
                'rank': i + 1,
            })

        # ---- RI-primary injection for branched alkanes ----
        try:
            from config import (RI_INJECT_ALKANES, RI_INJECT_TOL, RI_INJECT_MAX)
        except ImportError:
            RI_INJECT_ALKANES = False
        if RI_INJECT_ALKANES:
            results = self._inject_ri_alkanes(
                results, measured_ri, tol=RI_INJECT_TOL, max_inject=RI_INJECT_MAX)

        return results

    def search_peaks_raw(self, peaks, scan_list, top_n=5, verbose=True,
                         override_spectra=None):
        """Batch search + RI re-ranking (Goodner predictions included).

        override_spectra: optional list (per peak) of (mz_array, int_array) to
        search INSTEAD of that peak's raw apex spectrum — used to feed
        deconvolved 'pure' spectra of co-eluting peaks. None entries fall back
        to the raw apex scan."""
        all_matches = []
        for i, peak in enumerate(peaks):
            if override_spectra is not None and override_spectra[i] is not None:
                mz_arr, int_arr = override_spectra[i]
            else:
                scan = scan_list[peak['apex_idx']]
                mz_arr, int_arr = scan['mz'], scan['intensity']
            matches = self.search_raw_spectrum(
                mz_arr, int_arr, top_n, measured_ri=peak.get('ri_measured'))
            all_matches.append(matches)
            if verbose and (i + 1) % 50 == 0:
                print(f"    [NIST-RAW] {i+1}/{len(peaks)} peaks searched")
        return all_matches

    def search_peaks_binned(self, peaks, intensity_matrix, mz_bins,
                             top_n=5, verbose=True):
        """Batch search using binned spectra (clean sorted m/z, best NIST match)."""
        all_matches = []
        for i, peak in enumerate(peaks):
            spec = intensity_matrix[peak['apex_idx'], :]
            mask = spec > 0
            if not mask.any():
                all_matches.append([])
                continue
            mz_arr = mz_bins[mask]
            int_arr = spec[mask]
            matches = self.search_raw_spectrum(
                mz_arr, int_arr, top_n,
                measured_ri=peak.get('ri_measured'))
            all_matches.append(matches)
            if verbose and (i + 1) % 50 == 0:
                print(f"    [NIST-BINNED] {i+1}/{len(peaks)} peaks searched")
        return all_matches
