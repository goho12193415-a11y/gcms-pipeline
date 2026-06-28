"""
Step 6: Retention Time Alignment (Batch)
=========================================
Polynomial RT correction for multi-sample batches.
Uses greedy nearest-neighbor peak matching between samples.
"""
import numpy as np


def align_retention_times(
    sample_peaks: list,
    reference_peaks: list,
    rt_tolerance: float = 0.2,
    poly_degree: int = 1
) -> list:
    """Align retention times of sample peaks to a reference.

    Algorithm:
      1. Greedy nearest-neighbor matching of peaks within rt_tolerance
      2. Polynomial fit: sample_RT → reference_RT
      3. Apply RT transform to all sample peaks

    Args:
        sample_peaks: Peaks to align (list of peak dicts)
        reference_peaks: Reference peaks (same format)
        rt_tolerance: Initial matching window (min)
        poly_degree: Polynomial order (1=linear, 2=quadratic)

    Returns:
        Corrected sample peaks with added 'apex_rt_original' key
    """
    # 1. Match peaks (greedy nearest-neighbor)
    matched_pairs = []
    ref_used = set()

    for sp in sorted(sample_peaks, key=lambda x: x['apex_rt']):
        best_ref_idx = None
        best_dist = rt_tolerance
        for i, rp in enumerate(reference_peaks):
            if i in ref_used:
                continue
            dist = abs(sp['apex_rt'] - rp['apex_rt'])
            if dist < best_dist:
                best_dist = dist
                best_ref_idx = i
        if best_ref_idx is not None:
            matched_pairs.append(
                (sp['apex_rt'], reference_peaks[best_ref_idx]['apex_rt'])
            )
            ref_used.add(best_ref_idx)

    # Need at least poly_degree+1 points for fitting
    if len(matched_pairs) < poly_degree + 1:
        return sample_peaks  # Not enough points, return uncorrected

    # 2. Fit polynomial RT mapping
    sample_rts = np.array([p[0] for p in matched_pairs])
    ref_rts = np.array([p[1] for p in matched_pairs])
    coeffs = np.polyfit(sample_rts, ref_rts, poly_degree)
    rt_transform = np.poly1d(coeffs)

    # 3. Apply correction
    corrected = []
    for p in sample_peaks:
        p_copy = dict(p)
        p_copy['apex_rt_original'] = p['apex_rt']
        p_copy['apex_rt'] = float(rt_transform(p['apex_rt']))
        if 'start_rt' in p_copy:
            p_copy['start_rt_original'] = p_copy.get('start_rt')
            p_copy['start_rt'] = float(rt_transform(p_copy['start_rt_original']))
        if 'end_rt' in p_copy:
            p_copy['end_rt_original'] = p_copy.get('end_rt')
            p_copy['end_rt'] = float(rt_transform(p_copy['end_rt_original']))
        corrected.append(p_copy)

    return corrected
