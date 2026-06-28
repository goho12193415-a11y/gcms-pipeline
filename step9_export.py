"""
Step 9: Results Export v2 — Review-optimized Excel
====================================================
  - Top-3 candidates per peak (not just top-1)
  - Auto-confirm rules (RMF>900 + RI<20 + common food volatile)
  - Contaminant auto-tagging (siloxanes, phthalates, column bleed)
  - Area descending sort (biggest peaks first)
  - Color-coded: GREEN=auto-confirmed, YELLOW=review needed,
                 RED=contaminant, GRAY=low confidence
"""
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side


# ---- Contaminant patterns ----
CONTAMINANT_PATTERNS = [
    'silane', 'siloxane', 'dimethylsilanediol', 'silanediol',
    'phthalate', 'phthalic',
    'butylated hydroxytoluene', 'bht',
    'column bleed', 'cyclosiloxane',
    '2,4-di-tert-butylphenol', '2,6-di-tert-butyl',
    'diethyl phthalate', 'dibutyl phthalate', 'diisobutyl phthalate',
    'hexadecanoic acid', 'octadecanoic acid',  # common fatty acids (often contaminants in SPME)
]

# ---- Common food volatiles (auto-confirm candidates) ----
COMMON_FOOD_VOLATILES = [
    'hexanal', 'heptanal', 'octanal', 'nonanal', 'decanal',
    'benzaldehyde', 'phenylacetaldehyde',
    '2-heptenal', '2-octenal', '2-nonenal', '2-decenal',
    '2,4-heptadienal', '2,4-nonadienal', '2,4-decadienal',
    '2-pentylfuran', '2-ethylfuran',
    '1-octen-3-ol', '1-octen-3-one', '1-penten-3-ol',
    '2-heptanone', '2-octanone', '2-nonanone', '2-undecanone',
    'acetophenone', 'beta-ionone', '6,10,14-trimethylpentadecan-2-one',
    'limonene', 'alpha-pinene', 'beta-pinene', 'beta-caryophyllene',
    'linalool', 'alpha-terpineol', 'geraniol',
    '2-methylpyrazine', '2,5-dimethylpyrazine', '2,3,5-trimethylpyrazine',
    '2-ethylpyrazine', '2-ethyl-3,5-dimethylpyrazine',
    'dimethyl disulfide', 'dimethyl trisulfide', 'methional',
    'benzothiazole', 'indole',
    'phenol', 'guaiacol', '4-ethylguaiacol', 'eugenol',
    'hexanoic acid', 'octanoic acid', 'butanoic acid',
    'ethyl hexanoate', 'ethyl octanoate', 'ethyl decanoate',
    'gamma-butyrolactone', 'gamma-hexalactone',
    '2-ethyl-1-hexanol', '1-hexanol', '1-octanol', '1-nonanol',
    '2,2,4-trimethyl-1,3-pentanediol monoisobutyrate',
    '2,2,4-trimethyl-1,3-pentanediol diisobutyrate',
    'furan, 2-pentyl-', 'furan, 2-ethyl-',
    '2,4,6-trimethylpyridine', 'pyridine, 2,4,6-trimethyl-',
    'cyclocitral', 'beta-cyclocitral', 'photocitral',
    'farnesyl acetone', '(e,e)-farnesyl acetone',
    'alpha-isomethylionone', 'isomethylionone',
    '1-hexadecene', '1-octadecene',
    'dodecane', 'tridecane', 'tetradecane', 'pentadecane',
    'hexadecane', 'heptadecane', 'octadecane',
    '3,5-di-tert-butylphenol',
    '2-(2-butoxyethoxy)ethanol', 'ethanol, 1-(2-butoxyethoxy)-',
]


def _is_contaminant(name):
    """Check if compound name matches contaminant patterns."""
    n = name.lower()
    for pat in CONTAMINANT_PATTERNS:
        if pat in n:
            return True
    return False


def _is_common_food(name):
    """Check if compound is a common food volatile."""
    n = name.lower().strip()
    if n in COMMON_FOOD_VOLATILES:
        return True
    # Fuzzy match: if name contains a common food volatile as substring
    for fv in COMMON_FOOD_VOLATILES:
        if len(fv) > 10 and fv in n:
            return True
    return False


def _auto_confirm_status(matches, ri_tolerance=30):
    """
    Determine auto-review status for a peak.

    Returns: ('CONFIRMED', 'REVIEW', 'CONTAMINANT', 'LOW')
    """
    if not matches:
        return ('LOW', 'No match found')

    best = matches[0]
    name = best.get('name', '')
    rmf = best.get('rmf', 0)
    fmf = best.get('fmf', 0)
    ri_diff = best.get('ri_diff')

    # Contaminant check (highest priority)
    if _is_contaminant(name):
        return ('CONTAMINANT', 'Known contaminant / column bleed')

    coelution = best.get('coelution_flag', False)
    ri_ok = ri_diff is not None and ri_diff <= ri_tolerance

    # Auto-confirm: RMF>900 AND FMF>700 (both directions must agree, no co-elution)
    # AND (RI match < 30 OR common food volatile with no RI data)
    if rmf >= 900 and fmf >= 700 and not coelution:
        if ri_ok or (ri_diff is None and _is_common_food(name)):
            return ('CONFIRMED',
                    f'RMF={rmf:.0f} FMF={fmf:.0f} RI={"OK" if ri_ok else "N/A"} food')

    # Low confidence: co-elution or very unbalanced FMF/RMF
    if coelution or (fmf > 0 and rmf / max(fmf, 1) > 1.8):
        return ('LOW', f'Co-elution: RMF/FMF={rmf/max(fmf,1):.1f}')

    # Low confidence: RMF too low
    if rmf < 750 or fmf < 400:
        return ('LOW', f'Low: RMF={rmf:.0f} FMF={fmf:.0f}')

    # Needs review (middle ground)
    return ('REVIEW', f'Review: RMF={rmf:.0f} FMF={fmf:.0f}')


def compile_results(sample_name, integrated_peaks, identification_results,
                    quantification_results, metadata=None):
    """Compile results with top-3 candidates per peak.

    identification_results[i] = [match1, match2, match3, ...]
    """
    rows = []

    for i, peak in enumerate(integrated_peaks):
        matches = identification_results[i] if i < len(identification_results) else []
        quant = quantification_results[i] if i < len(quantification_results) else {}

        status, status_reason = _auto_confirm_status(matches)
        top = matches[0] if matches else {}

        # Peak info (shared across all candidates for this peak)
        base = {
            'Sample': sample_name,
            'Peak_No': i + 1,
            'RT_min': round(peak.get('apex_rt', 0), 3),
            'RT_start': round(peak.get('start_rt', 0), 3),
            'RT_end': round(peak.get('end_rt', 0), 3),
            'Area': round(peak.get('area', 0), 1),
            'Height': round(peak.get('height', 0), 1),
            'SN': round(peak.get('sn', 0), 1),
            'Percentage': round(quant.get('percentage', 0), 3) if quant.get('percentage') else '',
            'Status': status,
            'Status_Reason': status_reason,
            'User_Confirm': '',  # Empty column for user to fill
        }

        # Check if any known food volatile appears in top-5 (even if not rank 1)
        food_in_top5 = False
        for m in matches[:5]:
            if _is_common_food(m.get('name', '')):
                food_in_top5 = True
                break

        # Top-3 candidates as sub-rows
        for j, m in enumerate(matches[:3]):
            row = dict(base)
            row['Rank'] = j + 1
            row['Compound_Name'] = m.get('name', 'Unknown')
            row['CAS'] = m.get('cas', '')
            row['FMF'] = m.get('fmf', 0)
            row['RMF'] = m.get('rmf', 0)
            row['RI_WAX'] = m.get('ri_wax', '') if m.get('ri_wax') is not None else ''
            row['RI_DB5'] = m.get('ri_db5', '') if m.get('ri_db5') is not None else ''
            row['RI_Diff'] = m.get('ri_diff', '') if m.get('ri_diff') is not None else ''
            row['Coelution'] = 'Y' if m.get('coelution_flag') else ''
            row['Food_In_Top5'] = 'Y' if food_in_top5 else ''
            rows.append(row)

        # If no matches, still output one row
        if not matches:
            row = dict(base)
            row['Rank'] = ''
            row['Compound_Name'] = 'Unknown'
            row['CAS'] = ''
            row['FMF'] = ''
            row['RMF'] = ''
            row['RI_Diff'] = ''
            row['RI_Library'] = ''
            row['Coelution'] = ''
            rows.append(row)

    df = pd.DataFrame(rows)

    # Sort by Area descending (biggest peaks first, most reliable)
    if 'Area' in df.columns and not df.empty:
        df['_sort_area'] = pd.to_numeric(df['Area'], errors='coerce').fillna(0)
        df['_sort_rank'] = df.groupby('Peak_No')['_sort_area'].transform('max')
        # Stable sort: keep Peak_No ascending within same area
        df = df.sort_values(['_sort_rank', 'Peak_No', 'Rank'],
                            ascending=[False, True, True])
        df = df.drop(columns=['_sort_area', '_sort_rank'])

    return df


def export_to_excel(results_df, output_path, calibration_data=None):
    """Export to color-coded Excel with review workflow."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
        results_df.to_excel(writer, sheet_name='Results', index=False)

    # ---- Apply formatting with openpyxl ----
    from openpyxl import load_workbook
    wb = load_workbook(output_path)
    ws = wb['Results']

    # Styles
    green_fill = PatternFill(start_color='C6EFCE', end_color='C6EFCE', fill_type='solid')
    yellow_fill = PatternFill(start_color='FFEB9C', end_color='FFEB9C', fill_type='solid')
    red_fill = PatternFill(start_color='FFC7CE', end_color='FFC7CE', fill_type='solid')
    gray_fill = PatternFill(start_color='D9D9D9', end_color='D9D9D9', fill_type='solid')

    hdr_fill = PatternFill(start_color='2F5496', end_color='2F5496', fill_type='solid')
    hdr_font = Font(name='Calibri', size=10, bold=True, color='FFFFFF')
    body_font = Font(name='Calibri', size=10)
    thin_border = Border(
        left=Side(style='thin'), right=Side(style='thin'),
        top=Side(style='thin'), bottom=Side(style='thin')
    )
    center_align = Alignment(horizontal='center', vertical='center', wrap_text=True)

    # Freeze header
    ws.freeze_panes = 'A2'

    # Style header
    for col in range(1, ws.max_column + 1):
        cell = ws.cell(row=1, column=col)
        cell.font = hdr_font
        cell.fill = hdr_fill
        cell.alignment = center_align
        cell.border = thin_border

    # Column widths
    widths = {'A': 10, 'B': 6, 'C': 8, 'D': 8, 'E': 8, 'F': 12, 'G': 10, 'H': 6,
              'I': 8, 'J': 10, 'K': 14, 'L': 12, 'M': 30, 'N': 16, 'O': 6, 'P': 6,
              'Q': 8, 'R': 10, 'S': 6, 'T': 12}
    for col_letter, w in widths.items():
        if col_letter in [chr(65+i) for i in range(ws.max_column)]:
            ws.column_dimensions[col_letter].width = w

    # Find status column
    status_col = None
    for col in range(1, ws.max_column + 1):
        if ws.cell(row=1, column=col).value == 'Status':
            status_col = col
            break

    # Apply row colors based on Status
    for row in range(2, ws.max_row + 1):
        status = ws.cell(row=row, column=status_col).value if status_col else ''
        fill = gray_fill
        if status == 'CONFIRMED':
            fill = green_fill
        elif status == 'REVIEW':
            fill = yellow_fill
        elif status == 'CONTAMINANT':
            fill = red_fill

        for col in range(1, ws.max_column + 1):
            cell = ws.cell(row=row, column=col)
            cell.fill = fill
            cell.font = body_font
            cell.border = thin_border
            cell.alignment = center_align

        # Bold the Rank=1 row for each peak
        rank_cell = None
        for col in range(1, ws.max_column + 1):
            if ws.cell(row=1, column=col).value == 'Rank':
                rank_cell = ws.cell(row=row, column=col)
                break
        if rank_cell and rank_cell.value == 1:
            for col in range(1, ws.max_column + 1):
                ws.cell(row=row, column=col).font = Font(name='Calibri', size=10, bold=True)

    wb.save(output_path)

    # ---- Summary sheet ----
    n_total = results_df['Peak_No'].nunique()
    status_counts = results_df[results_df['Rank'].isin([1, '', None])]['Status'].value_counts()
    n_confirmed = status_counts.get('CONFIRMED', 0)
    n_review = status_counts.get('REVIEW', 0)
    n_contam = status_counts.get('CONTAMINANT', 0)
    n_low = status_counts.get('LOW', 0)

    print(f"  [Export] {output_path}")
    print(f"           {n_total} peaks: "
          f"GREEN={n_confirmed} YELLOW={n_review} RED={n_contam} GRAY={n_low}")
    print(f"           User reviews ~{n_review} peaks (skip {n_confirmed+n_contam+n_low} auto-classified)")
