"""
Step 7: Library Search v2 — Food sub-library + RI boost + RMF/FMF penalty
===========================================================================
Improvements over v1:
  1. Food volatile sub-library (65K vs 240K, 4x faster)
  2. RI match boosts combined score (was neutral)
  3. RMF >> FMF penalty (co-elution detection)
  4. Top-5 ion pre-search (was top-12)
"""
import json, pickle, re
import numpy as np
from collections import defaultdict
from pathlib import Path
from config import ALKANE_RTS_DB_WAX, RI_COLUMN


# ---- RI Calculation ----
def calc_ri_from_rt(rt, alkane_rts=None):
    if alkane_rts is None: alkane_rts = ALKANE_RTS_DB_WAX
    if not alkane_rts or rt < alkane_rts[0][0]:
        if len(alkane_rts) >= 2:
            a, b = alkane_rts[0], alkane_rts[1]
            return round(a[1] + (b[1]-a[1])/(b[0]-a[0]) * (rt-a[0]))
        return None
    if rt > alkane_rts[-1][0]:
        if len(alkane_rts) >= 2:
            a, b = alkane_rts[-2], alkane_rts[-1]
            return round(b[1] + (b[1]-a[1])/(b[0]-a[0]) * (rt-b[0]))
        return None
    for i in range(len(alkane_rts)-1):
        rt_n, ri_n = alkane_rts[i]; rt_n1, ri_n1 = alkane_rts[i+1]
        if rt_n <= rt <= rt_n1:
            return round((ri_n//100)*100 + 100*(rt-rt_n)/(rt_n1-rt_n))
    return None


# ---- Spectral Matching ----
def _weight_spectrum(int_arr, mz_bins):
    """Old pipeline formula: w = mz * sqrt(I / max_I) — better for small molecules."""
    w = np.zeros_like(int_arr, dtype=np.float64)
    max_i = int_arr.max()
    if max_i <= 0: return w
    mask = int_arr > 0
    w[mask] = mz_bins[mask] * np.sqrt(int_arr[mask] / max_i)
    return w

def compute_match_factor(sv, lv, mz_bins):
    """NIST MF/RMF with old pipeline weighting (n=1, m=0.5)."""
    ws = _weight_spectrum(sv, mz_bins)
    wl = _weight_spectrum(lv, mz_bins)
    ns, nl = np.dot(ws, ws), np.dot(wl, wl)
    if ns == 0 or nl == 0: return {'fmf': 0, 'rmf': 0}

    # Forward MF: all shared peaks
    dot = np.dot(ws, wl)
    fmf = int(round(999 * dot / np.sqrt(ns * nl)))
    fmf = max(0, min(999, fmf))

    # Reverse MF: only library peaks
    lib_mask = lv > 0
    ws_r = ws * lib_mask
    ns_r = np.dot(ws_r, ws_r)
    if ns_r == 0:
        rmf = 0
    else:
        shared = np.dot(ws_r, wl)
        rmf = int(round(999 * shared / np.sqrt(ns_r * nl)))
        rmf = max(0, min(999, rmf))

    return {'fmf': fmf, 'rmf': rmf}


# ---- SpectralLibrary with food sub-library support ----
class SpectralLibrary:
    def __init__(self):
        self.entries = []; self.bp_index = defaultdict(list)

    def load_json(self, library_path, mz_bins, min_ions=5):
        cache_path = str(Path(library_path).with_suffix('.pkl'))
        if Path(cache_path).exists():
            import os
            if os.path.getmtime(cache_path) > os.path.getmtime(library_path):
                print(f"  [LIB] Loading from cache: {Path(cache_path).name}")
                self._load_cache(cache_path, mz_bins)
                return

        print(f"  [LIB] Loading: {Path(library_path).name}")
        self._load_from_json(library_path, mz_bins, min_ions)

    def _load_cache(self, cache_path, mz_bins):
        with open(cache_path, 'rb') as f: c = pickle.load(f)
        n = len(c['names'])
        for i in range(n):
            self.entries.append({
                'name': c['names'][i], 'cas': c['cases'][i], 'formula': '',
                'ri': int(c['ris'][i]) if c['ris'][i] >= 0 else None,
                'spectrum': c['spectra'][i].astype(np.float32),
                'bp_idx': int(c['bp_indices'][i]),
                'ion_set': set(c['ion_sets'][i].tolist() if hasattr(c['ion_sets'][i], 'tolist') else list(c['ion_sets'][i])),
            })
        self.bp_index = defaultdict(list, {int(k): list(v) for k, v in c['bp_index'].items()})
        print(f"  [LIB] {len(self.entries):,} compounds ({len(self.bp_index):,} BP keys)")

    def _load_from_json(self, library_path, mz_bins, min_ions):
        with open(library_path, 'r', encoding='utf-8') as f: raw = json.load(f)
        n_mz = len(mz_bins)
        for entry in raw:
            if len(entry.get('mz_list', [])) < min_ions: continue
            vec = np.zeros(n_mz, dtype=np.float32)
            for mz, i in zip(entry['mz_list'], entry['intensity_list']):
                idx = np.argmin(np.abs(mz_bins - mz))
                if np.abs(mz_bins[idx] - mz) <= 0.5: vec[idx] += i
            if vec.max() > 0: vec = vec / vec.max() * 999.0
            top5 = np.argsort(vec)[-5:]
            ion_set = set(top5[vec[top5] > 0])
            bp = int(np.argmax(vec)) if vec.max() > 0 else 0
            self.entries.append({
                'name': entry.get('name', ''), 'cas': entry.get('cas', ''),
                'formula': entry.get('formula', ''), 'ri': entry.get('ri'),
                'spectrum': vec, 'bp_idx': bp, 'ion_set': ion_set,
            })
        for i, e in enumerate(self.entries):
            bp = e['bp_idx']
            for d in (-1, 0, 1):
                if 0 <= bp + d < n_mz: self.bp_index[bp + d].append(i)
        print(f"  [LIB] Indexed {len(self.entries):,} compounds ({len(self.bp_index):,} BP keys)")

    def pre_search(self, query_vec, max_candidates=300, min_shared=3):
        n_mz = len(query_vec)
        q_top5 = np.argsort(query_vec)[-5:]
        q_bp = int(np.argmax(query_vec))
        candidates = set()
        for d in (-1, 0, 1):
            if 0 <= q_bp + d < n_mz: candidates.update(self.bp_index.get(q_bp + d, []))
        if not candidates: return []
        q_ion_set = set()
        for idx in q_top5:
            if query_vec[idx] > 0: q_ion_set.update({idx-1, idx, idx+1})
        scored = []
        for idx in candidates:
            shared = sum(1 for i in self.entries[idx]['ion_set'] if i in q_ion_set)
            if shared >= min_shared: scored.append((shared, idx))
        scored.sort(reverse=True)
        return [idx for _, idx in scored[:max_candidates]]


def search_library(query_spectrum, lib, mz_bins, top_n=10,
                   max_candidates=500, min_rmf=700, ri_tolerance=30,
                   query_ri=None, ri_weight=0.05, coelution_penalty=True):
    """Search with RI boost + RMF>>FMF penalty."""
    candidates = lib.pre_search(query_spectrum, max_candidates=max_candidates)
    if not candidates: return []

    results = []
    for idx in candidates:
        e = lib.entries[idx]
        mf = compute_match_factor(query_spectrum, e['spectrum'], mz_bins)

        if mf['rmf'] < min_rmf: continue

        # RI scoring — strongest orthogonal evidence
        ri_diff = None
        ri_bonus = 0
        if query_ri is not None and e['ri'] is not None:
            ri_diff = abs(query_ri - e['ri'])
            if ri_diff > ri_tolerance * 2: continue  # Hard reject
            ri_score = max(0.0, 1.0 - ri_diff / ri_tolerance)
            ri_bonus = ri_score * 999 * ri_weight * 3.0  # 3x boost for RI match

        # No RI data → neutral (0.5 score)
        elif query_ri is not None and e['ri'] is None:
            ri_bonus = 0.5 * 999 * ri_weight  # neutral

        # RMF>>FMF penalty (co-elution: strong reverse, weak forward = mixed spectrum)
        rmf_fmf_ratio = mf['rmf'] / max(mf['fmf'], 1)
        co_penalty = 0
        if coelution_penalty and rmf_fmf_ratio > 1.5:
            co_penalty = min(250, (rmf_fmf_ratio - 1.0) * 150)  # Stronger penalty

        combined = mf['rmf'] * (1 - ri_weight) + ri_bonus - co_penalty
        combined = max(0, combined)

        results.append({
            'name': e['name'], 'cas': e['cas'],
            'formula': e.get('formula', ''),
            'fmf': mf['fmf'], 'rmf': mf['rmf'],
            'ri': e['ri'],
            'ri_diff': round(ri_diff, 1) if ri_diff is not None else None,
            'combined_score': round(combined, 1),
            'coelution_flag': rmf_fmf_ratio > 1.5,
        })

    results.sort(key=lambda x: x['combined_score'], reverse=True)
    return results[:top_n]
