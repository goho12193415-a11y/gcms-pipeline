"""
Step 4: Spectral Deconvolution
===============================
AMDIS-style ion-clustering deconvolution:
  1. Build 2D intensity matrix (scans × mz-bins)
  2. For each peak region, compute EIC correlation matrix
  3. Cluster correlated ions → pure component spectra
"""
import numpy as np
from scipy.stats import pearsonr


def build_intensity_matrix(scan_list: list, rt: np.ndarray,
                            mz_min: float = 35, mz_max: float = 500,
                            mz_step: float = 1.0) -> tuple:
    """Convert scan list to 2D intensity matrix.

    Rows = scan points (time), Columns = m/z bins.

    Args:
        scan_list: List of scan dicts with 'mz' and 'intensity' arrays
        rt: Retention time array
        mz_min, mz_max, mz_step: m/z binning parameters

    Returns:
        (intensity_matrix: np.ndarray [n_scans, n_mz], mz_bins: np.ndarray)
    """
    mz_bins = np.arange(mz_min, mz_max + mz_step, mz_step)
    n_scans = len(scan_list)
    n_mz = len(mz_bins)

    matrix = np.zeros((n_scans, n_mz), dtype=np.float32)

    for i, scan in enumerate(scan_list):
        if len(scan['mz']) == 0:
            continue
        # Assign each m/z to nearest bin
        bin_indices = np.round((scan['mz'] - mz_min) / mz_step).astype(int)
        valid = (bin_indices >= 0) & (bin_indices < n_mz)
        np.add.at(matrix[i], bin_indices[valid], scan['intensity'][valid])

    return matrix, mz_bins


def _cluster_correlated_ions(corr_matrix: np.ndarray, active_cols: np.ndarray,
                              threshold: float, min_ions: int) -> list:
    """Greedy clustering: group ions with pairwise correlation > threshold."""
    n = len(active_cols)
    used = [False] * n
    clusters = []

    for i in range(n):
        if used[i]:
            continue
        cluster = [i]
        used[i] = True
        for j in range(i + 1, n):
            if not used[j] and corr_matrix[i, j] >= threshold:
                cluster.append(j)
                used[j] = True
        if len(cluster) >= min_ions:
            clusters.append([active_cols[k] for k in cluster])

    return clusters if clusters else [list(active_cols)]


def deconvolve_peak_region(
    intensity_matrix: np.ndarray,
    mz_bins: np.ndarray,
    peak: dict,
    correlation_threshold: float = 0.80,
    min_ions: int = 3
) -> list:
    """Deconvolve a single chromatographic peak region.

    For each peak, extracts the time region and clusters ions
    whose EIC profiles are highly correlated (same compound).

    Args:
        intensity_matrix: Full 2D matrix [scans, mz_bins]
        mz_bins: m/z bin centers
        peak: Peak dict from Step 3 with start_idx, end_idx, apex_idx
        correlation_threshold: Pearson r threshold for ion clustering
        min_ions: Minimum ions per component

    Returns:
        List of deconvoluted components, each:
            {apex_idx, pure_spectrum: np.array, contributing_mz: np.array}
    """
    start = peak['start_idx']
    end = peak['end_idx'] + 1
    apex_rel = peak['apex_idx'] - start  # relative index within region

    region = intensity_matrix[start:end, :]  # (time_points, mz_bins)

    # Too short: return apex spectrum as-is
    if region.shape[0] < 3:
        spectrum = intensity_matrix[peak['apex_idx'], :].copy()
        active_mz = mz_bins[spectrum > 0]
        return [{
            'apex_idx': peak['apex_idx'],
            'pure_spectrum': spectrum,
            'contributing_mz': active_mz
        }]

    # Find active m/z columns (any signal in region)
    active_mask = region.max(axis=0) > 0
    active_cols = np.where(active_mask)[0]

    if len(active_cols) < min_ions:
        spectrum = intensity_matrix[peak['apex_idx'], :].copy()
        return [{
            'apex_idx': peak['apex_idx'],
            'pure_spectrum': spectrum,
            'contributing_mz': mz_bins[active_cols]
        }]

    eics = region[:, active_cols]  # (time, n_active)

    # Normalize each EIC (avoid intensity bias in correlation)
    eic_means = eics.mean(axis=0, keepdims=True)
    eic_stds = eics.std(axis=0, keepdims=True)
    eic_stds[eic_stds == 0] = 1.0
    eics_norm = (eics - eic_means) / eic_stds

    # Compute Pearson correlation matrix
    n_active = eics_norm.shape[1]
    corr_matrix = np.eye(n_active)
    for i in range(n_active):
        for j in range(i + 1, n_active):
            r, _ = pearsonr(eics_norm[:, i], eics_norm[:, j])
            if not np.isnan(r):
                corr_matrix[i, j] = corr_matrix[j, i] = r

    # Cluster correlated ions
    clusters = _cluster_correlated_ions(
        corr_matrix, active_cols, correlation_threshold, min_ions
    )

    # Reconstruct pure spectrum for each component
    result = []
    for comp_cols in clusters:
        apex_spectrum = np.zeros(len(mz_bins), dtype=np.float32)
        apex_time_idx = min(apex_rel, region.shape[0] - 1)
        for col in comp_cols:
            apex_spectrum[col] = region[apex_time_idx, col]

        result.append({
            'apex_idx': peak['apex_idx'],
            'pure_spectrum': apex_spectrum,
            'contributing_mz': mz_bins[comp_cols]
        })

    return result


def deconvolve_all_peaks(intensity_matrix: np.ndarray, mz_bins: np.ndarray,
                          peaks: list, **kwargs) -> list:
    """Run deconvolution on all detected peaks.

    Returns list-of-lists: deconv_results[peak_i] = [component, ...]
    """
    results = []
    for peak in peaks:
        components = deconvolve_peak_region(
            intensity_matrix, mz_bins, peak, **kwargs
        )
        results.append(components)
    return results
