"""
Step 3: Chromatographic Peak Detection
=======================================
TIC-based peak detection using derivative method (Genesis-equivalent).
Filters by S/N ratio and minimum peak width.
"""
import numpy as np
from scipy.signal import find_peaks, peak_widths
from scipy.ndimage import gaussian_filter1d


def _estimate_noise(signal: np.ndarray, quantile: float = 0.1) -> float:
    """Estimate noise: RMS of the lowest 10% intensity region.

    Equivalent to ICIS 'Area Noise Factor' estimation.
    """
    threshold = np.quantile(signal, quantile)
    noise_region = signal[signal <= threshold]
    if len(noise_region) < 3:
        return float(np.std(signal)) * 0.1
    return float(np.std(noise_region)) if np.std(noise_region) > 0 else 1.0


def detect_chromatographic_peaks(
    signal: np.ndarray,
    rt: np.ndarray,
    min_sn: float = 5.0,
    min_peak_width_scans: int = 3,
    smooth_sigma: float = 1.0,
    solvent_delay_min: float = 0.0
) -> list:
    """Detect chromatographic peaks in TIC signal.

    Algorithm:
      1. Gaussian pre-smoothing (σ = smooth_sigma scans)
      2. Derivative-based peak apex detection (scipy find_peaks)
      3. Peak width at half-height via peak_widths
      4. S/N filter using bottom-10% RMS noise estimate
      5. Solvent delay cut

    Args:
        signal: 1D TIC signal (preferably baseline-corrected)
        rt: Retention time array (same length as signal)
        min_sn: Minimum S/N ratio (default 5)
        min_peak_width_scans: Minimum peak width in scan points
        smooth_sigma: Gaussian sigma for pre-smoothing
        solvent_delay_min: Exclude peaks before this RT (min)

    Returns:
        List of peak dicts, each with:
            {apex_idx, apex_rt, apex_intensity,
             start_idx, end_idx, start_rt, end_rt,
             width_scans, sn}
    """
    # 1. Pre-smoothing
    smoothed = gaussian_filter1d(signal.astype(float), sigma=smooth_sigma)

    # 2. Peak detection (derivative-based)
    peak_indices, properties = find_peaks(
        smoothed,
        width=min_peak_width_scans,
        prominence=0  # Filter by S/N later
    )

    if len(peak_indices) == 0:
        return []

    # 3. Peak widths at half-height
    widths, width_heights, left_ips, right_ips = peak_widths(
        smoothed, peak_indices, rel_height=0.5
    )

    # 4. Noise estimate and S/N filter
    noise_level = _estimate_noise(smoothed)

    # 5. Build peak list
    peaks = []
    for i, apex_idx in enumerate(peak_indices):
        apex_intensity = float(smoothed[apex_idx])
        sn = apex_intensity / noise_level if noise_level > 0 else 0

        if sn < min_sn:
            continue

        start_idx = max(0, int(left_ips[i]))
        end_idx = min(len(signal) - 1, int(right_ips[i]))

        apex_rt = float(rt[apex_idx])
        if apex_rt < solvent_delay_min:
            continue

        peaks.append({
            'apex_idx': int(apex_idx),
            'apex_rt': apex_rt,
            'apex_intensity': apex_intensity,
            'start_idx': start_idx,
            'end_idx': end_idx,
            'start_rt': float(rt[start_idx]),
            'end_rt': float(rt[end_idx]),
            'width_scans': int(end_idx - start_idx),
            'sn': float(sn)
        })

    return sorted(peaks, key=lambda p: p['apex_rt'])
