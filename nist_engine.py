"""
NIST MS Search Engine v3 — Hybrid RI + MS Search
==================================================
For each peak:
  1. NIST MS search → top-10 candidates by spectral match
  2. RI database search → compounds within +/-50 RI of measured value
  3. Merge both lists — RI-matched compounds added even if MS rank is low
  4. Re-rank: compounds in BOTH lists get top priority
"""
import json, numpy as np
from pyms.Spectrum import MassSpectrum
import pyms_nist_search


class NISTSearchEngine:
    def __init__(self, lib_path=None, work_dir=None):
        if lib_path is None:
            lib_path = r"C:\Users\go ho\Desktop\MSSEARCH\mainlib"
        if work_dir is None:
            import tempfile, os
            work_dir = os.path.join(tempfile.gettempdir(), "nist_work")
            os.makedirs(work_dir, exist_ok=True)

        self.engine = pyms_nist_search.Engine(
            lib_path, pyms_nist_search.NISTMS_MAIN_LIB, work_dir,
        )
        self._load_ri_database()

    def _load_ri_database(self):
        """Load dual-column RI database."""
        self.ri_db5 = {}   # {name_lower: ri}
        self.ri_wax = {}   # {name_lower: ri}

        paths = [
            r"C:\Users\go ho\Desktop\gcms_pipeline_v2\library\ri_dual_column.json",
            r"C:\Users\go ho\Desktop\gcms_pipeline_v2\library\ri_nist_full.json",
            r"C:\Users\go ho\Desktop\gcms_pipeline_v2\library\ri_enriched.json",
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

    def _search_by_ri(self, measured_ri, tolerance=50):
        """Find all compounds in RI database within tolerance of measured RI.

        Returns list of {name, ri, delta}
        """
        if measured_ri is None:
            return []

        results = []
        for bucket in range((measured_ri - tolerance) // 10 * 10,
                            measured_ri + tolerance + 10, 10):
            for item in self.ri_lookup.get(bucket, []):
                name, ri = item[0], item[1]
                delta = abs(ri - measured_ri)
                if delta <= tolerance:
                    results.append({'name': name, 'ri': ri, 'delta': delta})

        # Sort by delta, deduplicate by name
        results.sort(key=lambda x: x['delta'])
        seen = set()
        unique = []
        for r in results:
            key = r['name'].lower()
            if key not in seen:
                seen.add(key)
                unique.append(r)
        return unique

    def search_raw_spectrum(self, mz_array, int_array, top_n=5,
                              measured_ri=None):
        """Search using raw high-precision m/z (no binning).

        Args:
            mz_array: numpy array of raw m/z values (0.1 Da precision)
            int_array: numpy array of raw intensities
        """
        if len(mz_array) == 0:
            return []

        # Normalize to base peak = 999
        max_i = int_array.max()
        if max_i <= 0:
            return []
        int_norm = (int_array / max_i * 999).astype(int)

        # Keep only top-200 most intense ions (NIST search limit)
        if len(mz_array) > 200:
            top_idx = np.argsort(int_norm)[-200:]
            mz_list = mz_array[top_idx].tolist()
            int_list = int_norm[top_idx].tolist()
        else:
            mz_list = mz_array.tolist()
            int_list = int_norm.tolist()

        # NIST search
        try:
            ms = MassSpectrum(mz_list, int_list)
            hits = self.engine.full_search_with_ref_data(ms)
        except Exception as e:
            return []

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
                'source': 'MS',
                'combined_score': rmf,
            })

        return results[:top_n]

    def search_spectrum(self, spectrum, mz_bins, top_n=5,
                         measured_ri=None):
        """NIST MS search (binned spectrum — legacy path)."""
        mask = spectrum > 0
        if not mask.any():
            return []

        mz_list = mz_bins[mask].tolist()
        max_i = spectrum[mask].max()
        int_list = (spectrum[mask] / max_i * 999).astype(int).tolist()

        # ---- Pass 1: NIST MS search ----
        try:
            ms = MassSpectrum(mz_list, int_list)
            hits = self.engine.full_search_with_ref_data(ms)
        except Exception as e:
            hits = []

        ms_results = []
        for i, hit in enumerate(hits[:top_n]):
            if not isinstance(hit, tuple) or len(hit) < 2:
                continue
            sr, ref = hit[0], hit[1]
            name = getattr(ref, 'name', '') or ''
            cas = getattr(ref, 'cas', '') or ''
            mf = int(getattr(sr, 'match_factor', 0) or 0)
            rmf = int(getattr(sr, 'reverse_match_factor', 0) or 0)

            # Look up RI: CAS first, then name
            ri_lit = None
            ri_delta = None
            if cas and cas in self.ri_by_cas:
                ri_lit = self.ri_by_cas[cas].get('DB-WAX')
            if ri_lit is None:
                ri_lit = self.ri_by_name.get(name.lower())
            if ri_lit and measured_ri:
                ri_delta = abs(measured_ri - ri_lit)

            ms_results.append({
                'name': name, 'cas': cas,
                'fmf': mf, 'rmf': rmf,
                'ri': ri_lit, 'ri_diff': ri_delta,
                'source': 'MS',
                'combined_score': rmf,
            })

        return ms_results[:top_n]

    def search_peaks_raw(self, peaks, scan_list, top_n=5, verbose=True):
        """Batch search using raw high-precision spectra (no binning)."""
        all_matches = []
        for i, peak in enumerate(peaks):
            scan = scan_list[peak['apex_idx']]
            matches = self.search_raw_spectrum(
                scan['mz'], scan['intensity'], top_n,
                measured_ri=peak.get('ri_measured'))
            all_matches.append(matches)
            if verbose and (i + 1) % 50 == 0:
                print(f"    [NIST-RAW] {i+1}/{len(peaks)} peaks searched")
        return all_matches

    def search_peaks_batch(self, peaks, intensity_matrix, mz_bins,
                            top_n=10, verbose=True):
        all_matches = []
        for i, peak in enumerate(peaks):
            if peak.get('enhanced_full') is not None:
                spec = peak['enhanced_full']
            else:
                spec = intensity_matrix[peak['apex_idx'], :]

            matches = self.search_spectrum(spec, mz_bins, top_n,
                                            measured_ri=peak.get('ri_measured'))
            all_matches.append(matches)

            if verbose and (i + 1) % 50 == 0:
                print(f"    [NIST+RI] {i+1}/{len(peaks)} peaks searched")

        return all_matches
