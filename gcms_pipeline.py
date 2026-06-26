#!/usr/bin/env python3
"""
GC-MS Auto-Identification Pipeline v1.0
========================================
Automated peak detection, spectral matching, and compound identification for
GC-MS data. Outputs formatted Excel reports with compound names, CAS numbers,
SI match factors, and integrated peak areas.

Usage:
    python gcms_pipeline.py sample.txt output.xlsx
    python gcms_pipeline.py --batch "data/*.txt" ./results/

Requirements: numpy, scipy, openpyxl

The bundled food_volatiles.msp library contains 29 common volatile compounds
with NIST-verified reference EI-MS spectra. Additional spectra can be appended
in MSP format. For extended coverage, optionally load AMDIS/NIST MSL libraries.

Author: GC-MS Pipeline Contributors
License: MIT
"""

import os, sys, re, json, argparse
import numpy as np
from scipy.signal import savgol_filter, find_peaks
from scipy.ndimage import minimum_filter1d, gaussian_filter1d
from scipy.integrate import simpson

# ---- Config ----
AMDIS_LIB_DIR = r"C:\NIST26-EI-DEMO\AMDIS32\LIB"
SI_THRESHOLD = 650

class Config:
    threshold = 650
PEAK_DISTANCE = 10
PROMINENCE_FACTOR = 0.3

# Background markers to filter
AIR_IONS = {28, 32, 40, 44, 18, 17}
BLEED_IONS = {73, 147, 207, 267, 281, 355}
PHTHALATE_IONS = {149, 167, 279}
# Known contaminants (not sample-derived)
CONTAMINANTS = {'2,4-Di-tert-butylphenol', 'Butylated hydroxytoluene',
                'Diethyl phthalate', 'Dibutyl phthalate'}

# RT sanity ranges (approximate, for DB-5 MS 50-min program at ~5 °C/min from 40°C)
# Compounds eluting outside their expected range are false positives
RT_RANGES = {
    '2-Heptanol': (8, 22),        # C7 alcohol, BP ~160°C
    '2-Methyl-1-butanol': (6.5, 15), # C5 alcohol, BP ~129°C, DB-5 RT > 6.5 min
    '2-Heptanone': (4, 22),        # C7 ketone
    '1-Heptanol': (6, 22),         # C7 alcohol
    '2-Ethylhexanol': (6, 24),     # C8 alcohol
    '2-Nonanol': (10, 28),         # C9 alcohol
}


# ============================================================
# MSP LIBRARY PARSER
# ============================================================
def parse_msp(filepath):
    """Parse MSP (Mass Spectral Peak) format library.
    Entries are delimited by 'Name:' headers.

    Example entry:
        Name: Hexanal
        Formula: C6H12O
        CASNO: 66-25-1
        Num Peaks: 18
        44 999; 56 450; 41 400; ...
    """
    compounds = []
    with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
        content = f.read()

    # Split on "Name:" to get individual entries
    # First entry starts with 'Name:', rest start with '\nName:'
    raw_entries = []
    idx = 0
    while True:
        npos = content.find('Name:', idx)
        if npos == -1:
            break
        next_npos = content.find('\nName:', npos + 1)
        if next_npos == -1:
            raw_entries.append(content[npos:].strip())
            break
        raw_entries.append(content[npos:next_npos].strip())
        idx = next_npos + 1

    for entry in raw_entries:
        lines = entry.split('\n')
        name = ''
        formula = ''
        cas = ''
        peaks = []

        for i, line in enumerate(lines):
            line = line.strip()
            if not line:
                continue

            if line.startswith('Name: '):
                name = line[6:].strip()
            elif line.startswith('Formula: '):
                formula = line[9:].strip()
            elif line.startswith('CASNO: '):
                cas = line[7:].strip()
            elif line.startswith('Num Peaks:'):
                continue
            else:
                # Peak data: "mz int; mz int; ..."
                for part in line.replace(';', ' ').split():
                    part = part.strip()
                    # This is a bit tricky - we need to read mz,int pairs
                    pass
                # Better approach: use regex to find numbers
                import re as _re
                numbers = _re.findall(r'(\d+)\s+(\d+)', line)
                for mz_str, int_str in numbers:
                    peaks.append((int(mz_str), int(int_str)))

        if name and len(peaks) >= 5:
            top8 = [m for m, i in sorted(peaks, key=lambda x: -x[1])[:8]]
            fuzzy = set(top8)
            for m in top8: fuzzy.update({m-1, m+1})
            compounds.append({
                'name': name, 'formula': formula, 'cas': cas, 'peaks': peaks,
                'top8': top8, 'top8_fuzzy': fuzzy
            })

    return compounds


# ============================================================
# MSL LIBRARY PARSER
# ============================================================
def parse_msl(filepath):
    """Parse AMDIS/NIST MSL format library."""
    compounds = []
    with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
        content = f.read()

    entries = content.split('\nNAME: ')
    for i, entry in enumerate(entries):
        if i == 0 and entry.startswith('NAME: '):
            entry = entry[6:]
        if not entry.strip():
            continue

        lines = entry.strip().split('\n')
        name = lines[0].strip()
        formula = ''
        cas = ''
        peaks = []

        for line in lines[1:]:
            line = line.strip()
            if line.startswith('FORM:'):
                formula = line[5:].strip()
            elif line.startswith('CASNO:'):
                cas = line[6:].strip()
            else:
                m = re.findall(r'\(\s*(\d+)\s+(\d+)\s*\)', line)
                for mz_str, int_str in m:
                    peaks.append((int(mz_str), int(int_str)))

        if name and len(peaks) >= 5:
            top8 = [m for m, i in sorted(peaks, key=lambda x: -x[1])[:8]]
            fuzzy = set(top8)
            for m in top8: fuzzy.update({m-1, m+1})
            compounds.append({
                'name': name, 'formula': formula, 'cas': cas, 'peaks': peaks,
                'top8': top8, 'top8_fuzzy': fuzzy
            })

    return compounds


# ============================================================
# LIBRARY LOADER
# ============================================================
def load_libraries(lib_dir=None):
    """Load spectral libraries: bundled MSP + optional AMDIS MSL."""
    all_compounds = []

    # Primary: NIST main library (if available), fallback to food_volatiles.msp
    here = os.path.dirname(os.path.realpath(__file__))
    nist_path = os.path.join(here, 'nist_mainlib.msp')
    fallback_path = os.path.join(here, 'food_volatiles.msp')

    if os.path.exists(nist_path):
        comps = parse_msp(nist_path)
        all_compounds.extend(comps)
        print(f"[LIB] nist_mainlib.msp: {len(comps)} compounds (NIST 2014)")
    elif os.path.exists(fallback_path):
        comps = parse_msp(fallback_path)
        all_compounds.extend(comps)
        print(f"[LIB] food_volatiles.msp: {len(comps)} compounds")

    # Secondary: AMDIS MSL (only if --with-msl flag is set)
    if lib_dir and os.path.isdir(lib_dir) and getattr(Config, 'use_msl', False):
        for fname in ['NISTFF.MSL', 'NISTEPA.MSL']:
            fp = os.path.join(lib_dir, fname)
            if os.path.isfile(fp):
                comps = parse_msl(fp)
                all_compounds.extend(comps)
                print(f"[LIB] {fname}: {len(comps)} compounds (MSL)")
    elif lib_dir and os.path.isdir(lib_dir):
        print("[LIB] MSL libraries available but not loaded (use --with-msl to enable)")

    # Build base-peak index for ultrafast pre-search (Thermo-style)
    bp_index = {}
    for i, comp in enumerate(all_compounds):
        bp = max(comp['peaks'], key=lambda x: x[1])[0]
        bp = int(round(bp))
        if bp not in bp_index: bp_index[bp] = []
        bp_index[bp].append(i)

    print(f"[LIB] Total: {len(all_compounds)} compounds, indexed by {len(bp_index)} base peaks")
    return all_compounds, bp_index


# ============================================================
# SPECTRAL MATCHING
# ============================================================
def spectral_similarity(observed, reference):
    """NIST mass-weighted cosine with contaminant ion pre-filtering.

    - Filters known contaminant m/z before matching (siloxanes, phthalates)
    - Weights ions by sqrt(intensity) × m/z^1.3 (NIST convention)
    - Forward + reverse search with base peak gating
    """
    # Filter contaminant ions from observed spectrum before matching
    CONTAMINANT_MZ = {73, 147, 207, 221, 267, 281, 295, 327, 341, 355, 429,
                      149, 167, 279, 57, 71, 85, 99, 113}  # Siloxane + alkane + phthalate
    obs = {}
    for m, i in observed:
        mz = int(round(m))
        if mz not in CONTAMINANT_MZ and mz > 25:
            obs[mz] = max(obs.get(mz, 0), int(i))

    ref = {}
    for m, i in reference:
        mz = int(round(m))
        ref[mz] = max(ref.get(mz, 0), int(i))

    if len(obs) < 5 or len(ref) < 5:
        return 0

    ref_bp_mz = max(ref, key=ref.get)
    obs_bp_mz = max(obs, key=obs.get)
    obs_top3 = sorted(obs.keys(), key=lambda x: -obs[x])[:3]

    # Gate: ref base peak must be in observed top 3 (±2 Da, relaxed for co-elution)
    if not any(abs(ref_bp_mz - o) <= 2 for o in obs_top3):
        return 0

    # Normalize to 0-1 range
    obp = max(obs.values()); onorm = {m: v/obp for m, v in obs.items()}
    rbp_val = max(ref.values()); rnorm = {m: v/rbp_val for m, v in ref.items()}

    # NIST weight: W = sqrt(I) × mz^1.3
    def nist_weight(mz, intensity):
        return np.sqrt(max(intensity, 0.001)) * (mz ** 1.3)

    # Shared ions (within 1 Da tolerance)
    rkeys = set(rnorm.keys())
    shared = set()
    for m in onorm:
        for d in (-1, 0, 1):
            if m+d in rkeys:
                shared.add(m)
                break

    if len(shared) < 6:
        return 0

    # Mass-weighted cosine
    fnum = 0.0; fden_o = 0.0; fden_r = 0.0
    for mz in shared:
        oi = onorm.get(mz, 0.001); ri = rnorm.get(mz, 0.001)
        wo = nist_weight(mz, oi); wr = nist_weight(mz, ri)
        fnum += wo * wr; fden_o += wo*wo; fden_r += wr*wr

    if fden_o == 0 or fden_r == 0: return 0
    score = fnum / (np.sqrt(fden_o) * np.sqrt(fden_r))

    # Reverse coverage: what fraction of REF ions found in OBS?
    rev_hits = sum(1 for rmz in rnorm if any(abs(rmz+d-omz)<=1 for omz in onorm for d in (-1,0,1)))
    rev_cov = rev_hits / max(len(rnorm), 1)
    if rev_cov < 0.4: return 0

    # BP bonus
    bp_bonus = 1.0 if abs(ref_bp_mz - obs_bp_mz) <= 1 else (0.5 if abs(ref_bp_mz - obs_bp_mz) <= 3 else 0)
    if bp_bonus == 0: return 0

    return min(999, int((score*0.50 + rev_cov*0.30 + bp_bonus*0.20) * 999))


# ============================================================
# DATA PARSER
# ============================================================
def parse_thermo_txt(filepath):
    """Parse Thermo Xcalibur GC-MS text export."""
    rts, tics, spectra = [], [], []
    rt, tic, pkts, active = None, None, [], False

    with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith('ScanHeader'):
                if active and rt is not None:
                    rts.append(rt)
                    tics.append(tic or 0)
                    spectra.append(pkts)
                rt, tic, pkts, active = None, None, [], True
                continue
            if not active:
                continue
            if rt is None:
                m = re.search(r'start_time = ([\d.]+)', s)
                if m:
                    rt = float(m.group(1))
            m = re.search(r'integ_intens = ([\d.]+)', s)
            if m:
                tic = float(m.group(1))
            m = re.search(
                r'Packet # \d+, intensity = ([\d.eE+-]+), mass/position = ([\d.]+)',
                s
            )
            if m:
                pkts.append((float(m.group(2)), float(m.group(1))))
        if active and rt is not None:
            rts.append(rt)
            tics.append(tic or 0)
            spectra.append(pkts)

    return np.array(rts), np.array(tics), spectra


# ============================================================
# PEAK DETECTION
# ============================================================
def detect_peaks(rts, tics, spectra):
    """Detect peaks, filter background, integrate areas, extract spectra."""
    # Smooth
    tics_s = savgol_filter(tics, 31, 2)
    # Baseline
    bl = minimum_filter1d(tics_s, 401)
    bl = gaussian_filter1d(bl.astype(float), 80)
    corr = tics_s - bl
    corr[corr < 0] = 0

    # Find peaks
    thresh = np.percentile(corr[corr > 0], 50)
    idxs, _ = find_peaks(corr, height=thresh, distance=PEAK_DISTANCE,
                          prominence=thresh * PROMINENCE_FACTOR)
    min_tic = np.median(tics) * 1.2
    idxs = [i for i in idxs if tics[i] > min_tic]

    results = []
    for idx in idxs:
        # Integration bounds
        left = idx
        while left > 0 and corr[left] > corr[left - 1]:
            left -= 1
        right = idx
        while right < len(corr) - 1 and corr[right] > corr[right + 1]:
            right += 1

        x = rts[left:right + 1]
        y = tics[left:right + 1]
        base = np.linspace(tics[left], tics[right], len(y))
        yc = y - base
        yc[yc < 0] = 0
        area = float(simpson(yc, x)) if len(x) >= 3 else 0

        # Spectrum at apex
        raw = spectra[idx]
        top = sorted(raw, key=lambda x: -x[1])[:40]

        results.append({
            'rt': round(rts[idx], 4),
            'tic': round(tics[idx], 0),
            'area': round(area, 0),
            'spectrum': top,
            'scan': int(idx),
        })

    return results


# ============================================================
# BACKGROUND FILTER
# ============================================================
def classify_peak(mz_int_pairs):
    """3-tier classification:
    L1: Natural products (aldehydes, ketones, alcohols, esters, terpenes, etc.)
    L2: Suspected contaminants (branched alkanes, siloxanes, glycol ethers, phthalates)
    L3: Confirmed background (air/water, column bleed)
    Returns (tier, label, keep_for_matching)
    """
    mz_list = [int(round(m)) for m, i in mz_int_pairs[:15]]
    bp = mz_list[0] if mz_list else 0

    # == L3: Confirmed background (ignore entirely) ==
    # Air/water: dominated by m/z 18, 28, 32, 40, 44
    air_count = sum(1 for m in AIR_IONS if m in mz_list[:8])
    if air_count >= 5:
        return ('L3', 'Air/Water', False)

    # Column bleed: siloxane pattern m/z 73 + 147 + 207
    bleed_count = sum(1 for m in BLEED_IONS if m in mz_list)
    if bleed_count >= 2:
        return ('L3', 'Column Bleed', False)

    # == L2: Suspected contaminants (keep but flag) ==
    # Phthalate plasticizers: m/z 149 dominant
    if 149 in mz_list[:5]:
        return ('L2', 'Plasticizer', True)

    # Branched alkanes: BP=57 with strong 43,71,85 pattern
    if bp == 57 and all(m in mz_list for m in [43, 71, 85]):
        return ('L2', 'Branched Alkane', True)

    # Siloxane derivatives: BP=73 without full column bleed
    if bp == 73:
        return ('L2', 'Siloxane Derivative', True)

    # Glycol ethers / PEG: BP=45 with large fragment gaps
    if bp == 45 and len(mz_list) < 12:
        return ('L2', 'Glycol/PEG Derivative', True)

    # == L1: Natural products ==
    return ('L1', 'Natural Product', True)


# ============================================================
# MAIN PROCESSING
# ============================================================
def process_sample(filepath, library, bp_index, output_path):
    """Run full pipeline on one GC-MS data file."""
    name = os.path.basename(filepath)
    print(f"\n[Processing] {name}")

    # Parse
    rts, tics, spectra = parse_thermo_txt(filepath)
    print(f"  Scans: {len(rts)}, RT: {rts[0]:.1f}-{rts[-1]:.1f} min")

    # Detect & integrate
    peaks = detect_peaks(rts, tics, spectra)
    print(f"  Raw peaks: {len(peaks)}")

    # Identify
    results = []
    for p in peaks:
        raw = p['spectrum']
        # Build clean spectrum (m/z > 25, positive intensity)
        clean = [(int(round(m)), int(i)) for m, i in raw
                  if int(round(m)) > 25 and i > 100]
        if len(clean) < 5:
            continue

        # Classify (3-tier system)
        tier, peak_label, keep_for_matching = classify_peak(clean[:15])
        if not keep_for_matching:
            continue  # Only skip L3 (air/bleed)

        # For L2 peaks, skip the strict air filter
        if tier == 'L1':
            # Remove air ions for cleaner matching
            air_count = sum(1 for m, i in clean[:12] if m in AIR_IONS)
            if air_count >= 6:
                continue
            clean2 = [(m, i) for m, i in clean if m not in AIR_IONS]
            if len(clean2) >= 6:
                clean = clean2

        # Match against library (with presearch + tier-based threshold)
        # Presearch: filter by top 10 most abundant ions first
        effective_threshold = 700 if tier == 'L1' else (Config.threshold if Config.threshold > 850 else 850)

        # Thermo-style indexed pre-search: only check compounds whose base peak matches
        # observed spectrum's top 3 most abundant ions (±1 Da)
        obs_bp = int(round(max(clean, key=lambda x: x[1])[0]))
        obs_top3 = [int(round(m)) for m, i in sorted(clean, key=lambda x: -x[1])[:3]]

        # Collect candidate indices from bp_index
        candidates = set()
        for bp_candidate in [obs_bp-1, obs_bp, obs_bp+1]:
            if bp_candidate in bp_index:
                candidates.update(bp_index[bp_candidate])

        matches = []
        for idx in candidates:
            comp = library[idx]
            si = spectral_similarity(clean, comp['peaks'])
            if si >= effective_threshold:
                matches.append((si, comp))

        if not matches:
            continue

        # Filter out matches that are impossible in aqueous HS-SPME samples
        # These compounds either hydrolyze or are too volatile to retain
        aqueous_impossible = {
            'Ethyl acetate', 'Hexyl acetate', 'Bornyl acetate',
            '1-Octen-3-ol, acetate', 'Octanoic acid, ethyl ester',
            'Diethyl phthalate', 'Dibutyl phthalate',
        }
        filtered_matches = []
        for si, comp in matches:
            if comp['name'] in aqueous_impossible and si < 950:
                continue  # Reject unless extremely high confidence
            filtered_matches.append((si, comp))

        if not filtered_matches:
            continue

        # Check for phthalate pattern in observed spectrum
        has_phthalate_149 = any(m == 149 for m, i in clean[:20])
        has_phthalate_57 = any(m == 57 for m, i in clean[:20])
        if has_phthalate_149 and has_phthalate_57:
            continue  # Phthalate contaminant, not a sample compound

        matches = filtered_matches
        matches.sort(key=lambda x: -x[0])
        si, comp = matches[0]

        # Filter known contaminants
        if comp['name'] in CONTAMINANTS:
            continue

        # RT sanity check
        rt_range = RT_RANGES.get(comp['name'])
        if rt_range and (p['rt'] < rt_range[0] or p['rt'] > rt_range[1]):
            continue  # Compound eluting at physically impossible time

        confidence = 'High' if si >= 850 else ('Medium' if si >= 750 else 'Low')

        results.append({
            'rt': p['rt'],
            'area': p['area'],
            'tic': p['tic'],
            'compound': comp['name'],
            'cas': comp.get('cas', ''),
            'formula': comp.get('formula', ''),
            'si': si,
            'confidence': confidence,
            'tier': tier,
            'peak_type': peak_label,
        })

    # Deduplicate
    results.sort(key=lambda x: -x['si'])
    seen = set()
    final = []
    for r in results:
        key = (round(r['rt'], 1), r['compound'])
        if key not in seen:
            seen.add(key)
            final.append(r)
    final.sort(key=lambda x: x['rt'])

    # Report
    _write_excel(final, output_path, name)

    n_high = sum(1 for r in final if r['confidence'] == 'High')
    n_med = sum(1 for r in final if r['confidence'] == 'Medium')
    print(f"  Identified: {len(final)} (High:{n_high} Med:{n_med})")
    return final


def _write_excel(results, path, sample_name):
    """Generate formatted Excel report."""
    import openpyxl
    from openpyxl.styles import Font, Alignment, PatternFill, Border, Side

    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = sample_name[:31]

    hf = Font(name='Arial', size=10, bold=True, color='FFFFFF')
    hfill = PatternFill('solid', fgColor='2F5496')
    gf = PatternFill('solid', fgColor='C6EFCE')
    yf = PatternFill('solid', fgColor='FFEB9C')
    bd = Border(left=Side('thin'), right=Side('thin'),
                top=Side('thin'), bottom=Side('thin'))
    df = Font(name='Arial', size=10)
    da = Alignment(horizontal='center', vertical='center')

    hdrs = ['No.', 'RT(min)', 'Compound', 'CAS', 'SI', 'Confidence', 'Tier',
            'Peak Type', 'Peak Area', 'Peak Height', 'Formula']
    widths = [5, 10, 35, 15, 7, 10, 6, 14, 16, 16, 16]

    for c, (h, w) in enumerate(zip(hdrs, widths), 1):
        cell = ws.cell(row=1, column=c, value=h)
        cell.font = hf; cell.fill = hfill
        cell.alignment = Alignment(horizontal='center', vertical='center', wrap_text=True)
        cell.border = bd
        ws.column_dimensions[openpyxl.utils.get_column_letter(c)].width = w

    l2_fill = PatternFill('solid', fgColor='FFC000')  # Orange for L2 contaminants
    l3_fill = PatternFill('solid', fgColor='FF9999')  # Red for L3 background

    for i, r in enumerate(results, 2):
        vals = [i - 1, r['rt'], r['compound'], r['cas'], r['si'],
                r['confidence'], r.get('tier', ''), r.get('peak_type', ''),
                r['area'], r['tic'], r.get('formula', '')]
        for c, v in enumerate(vals, 1):
            cell = ws.cell(row=i, column=c, value=v)
            cell.font = df; cell.alignment = da; cell.border = bd

        # Color coding: L1=green/yellow, L2=orange, L3=red
        tier = r.get('tier', 'L1')
        if tier == 'L2':
            fill = l2_fill
        elif tier == 'L3':
            fill = l3_fill
        else:
            fill = gf if r['confidence'] == 'High' else (yf if r['confidence'] == 'Medium' else None)
        if fill:
            for c in range(1, len(hdrs) + 1):
                ws.cell(row=i, column=c).fill = fill

    ws.auto_filter.ref = f'A1:{openpyxl.utils.get_column_letter(len(hdrs))}{len(results)+1}'
    ws.freeze_panes = 'A2'
    wb.save(path)


# ============================================================
# CLI
# ============================================================
def main():
    ap = argparse.ArgumentParser(description='GC-MS Auto-Identification Pipeline')
    ap.add_argument('input', help='Input TXT file or glob pattern with --batch')
    ap.add_argument('output', help='Output Excel path or directory with --batch')
    ap.add_argument('--lib', default=AMDIS_LIB_DIR, help='AMDIS LIB directory path')
    ap.add_argument('--batch', action='store_true', help='Batch mode')
    ap.add_argument('--with-msl', action='store_true', help='Load AMDIS MSL libraries')
    ap.add_argument('--threshold', type=int, default=650,
                    help='Minimum SI threshold (default: 650)')
    args = ap.parse_args()

    Config.threshold = args.threshold
    Config.use_msl = args.with_msl

    print("=" * 60)
    print("GC-MS Auto-Identification Pipeline v1.0")
    print("=" * 60)

    library, bp_index = load_libraries(args.lib)

    if args.batch:
        import glob
        files = glob.glob(args.input)
        os.makedirs(args.output, exist_ok=True)
        for f in sorted(files):
            base = os.path.splitext(os.path.basename(f))[0]
            out = os.path.join(args.output, f'{base}_identified.xlsx')
            process_sample(f, library, bp_index, out)
    else:
        process_sample(args.input, library, bp_index, args.output)

    print("\nDone.")


if __name__ == '__main__':
    main()
