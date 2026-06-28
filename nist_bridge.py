"""
NIST MS Search Bridge
=====================
Export peak spectra as MSP file for batch search in NIST MS Search,
then parse the search results back.

Workflow:
  1. Pipeline runs → generates `sample_spectra.msp`
  2. User opens in NIST MS Search → batch search → export results to text
  3. Pipeline reads the text → merges with integration data
"""
import numpy as np
from pathlib import Path


def export_peaks_to_msp(peaks, intensity_matrix, mz_bins, rt_array,
                         output_path, sample_name="sample",
                         ri_column="DB-WAX"):
    """Export all peak spectra as an MSP file for NIST MS Search.

    Format follows NIST MSP specification:
      Name: peak_RT_min
      RI: retention_index
      Num Peaks: N
      mz intensity; mz intensity; ...

    NIST MS Search will read this and batch-search all entries.
    """
    lines = []
    for i, peak in enumerate(peaks):
        # Get enhanced spectrum if available, else apex
        if peak.get('enhanced_full') is not None:
            spec = peak['enhanced_full']
        else:
            spec = intensity_matrix[peak['apex_idx'], :]

        # Extract non-zero m/z-intensity pairs
        mask = spec > 0
        mz_vals = mz_bins[mask]
        int_vals = spec[mask]

        # Normalize to base peak = 999
        if int_vals.max() > 0:
            int_vals = int_vals / int_vals.max() * 999

        # Sort by m/z ascending (NIST standard)
        order = np.argsort(mz_vals)
        peaks_str = '; '.join(
            f"{int(mz_vals[j])} {int(int_vals[j])}"
            for j in order if int_vals[j] >= 1
        )

        ri = peak.get('ri_measured', '')
        ri_str = f"\nRI: {ri}" if ri else ""

        lines.append(
            f"Name: Peak_{i+1}_RT_{peak['apex_rt']:.3f}"
            f"{ri_str}"
            f"\nNum Peaks: {len(order)}"
            f"\n{peaks_str}\n"
        )

    with open(output_path, 'w') as f:
        f.write('\n'.join(lines))

    print(f"  [MSP] Exported {len(peaks)} spectra to {output_path}")
    return output_path


def parse_nist_export(txt_path, peaks, rt_tolerance=0.1):
    """Parse NIST MS Search text export results.

    Expected format (NIST MS Search "Report" output):
      Hit list for Peak_X_RT_Y.YYY
      ...
      Rank  Name                    MF   RMF   Prob   CAS        RI
      ----  ---------------------  ----  ----  -----  ---------  ----
      1     Hexanal                 900   910   95.2   66-25-1    1083
      2     ...

    Returns: list of match lists (one per peak)
      [{rank, name, cas, mf, rmf, prob, ri}, ...]
    """
    with open(txt_path, 'r', encoding='utf-8', errors='replace') as f:
        content = f.read()

    # Split by "Hit list for" blocks
    import re
    blocks = re.split(r'Hit list for (.+)', content)

    results = []
    current_peak_name = None

    for i, block in enumerate(blocks):
        if i == 0 and not block.strip():
            continue
        if i % 2 == 1:  # Peak name line
            current_peak_name = block.strip()
        else:  # Hit list
            # Parse peak number from name: Peak_X_RT_Y.YYY
            match = re.match(r'Peak_(\d+)_RT_([\d.]+)', current_peak_name or '')
            if not match:
                continue

            peak_idx = int(match.group(1)) - 1
            peak_rt = float(match.group(2))

            # Find matching peak
            if peak_idx < len(peaks):
                peak = peaks[peak_idx]
            else:
                # Search by RT
                peak = None
                for p in peaks:
                    if abs(p['apex_rt'] - peak_rt) < rt_tolerance:
                        peak = p
                        break

            if peak is None:
                continue

            # Parse hit lines
            hits = []
            for line in block.split('\n'):
                line = line.strip()
                if not line or line.startswith('--') or line.startswith('Rank'):
                    continue
                # Parse: Rank  Name  MF  RMF  Prob  CAS  RI
                parts = line.split()
                if len(parts) < 3:
                    continue
                try:
                    rank = int(parts[0])
                except ValueError:
                    continue

                # Name is everything until we hit a number (MF)
                name_parts = []
                for p in parts[1:]:
                    try:
                        float(p)
                        break  # Hit a number = end of name
                    except ValueError:
                        name_parts.append(p)
                name = ' '.join(name_parts)

                # Remaining numbers: MF, RMF, Prob, CAS?, RI?
                remaining = parts[1 + len(name_parts):]
                mf = int(remaining[0]) if len(remaining) > 0 else 0
                rmf = int(remaining[1]) if len(remaining) > 1 else 0
                prob = float(remaining[2]) if len(remaining) > 2 else 0
                cas = remaining[3] if len(remaining) > 3 and '-' in remaining[3] else ''
                ri = int(remaining[-1]) if len(remaining) > 4 else None

                hits.append({
                    'rank': rank, 'name': name,
                    'cas': cas, 'fmf': mf, 'rmf': rmf,
                    'prob': prob, 'ri': ri
                })

            results.append({
                'peak_idx': peak_idx if peak_idx < len(peaks) else None,
                'peak_rt': peak_rt,
                'hits': hits
            })

    return results
