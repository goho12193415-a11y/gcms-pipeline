"""
Layer 4: Spectrum Enhancement (Xcalibur's key accuracy secret)
================================================================
Per GCMS_Complete_Workflow_v2.md:
  - Multi-scan averaging at peak apex (3-5 scans)
  - Subtract background spectrum from peak flanks
  - SNR improvement of √n, removes column bleed / matrix ions
"""
import numpy as np


def enhance_spectrum(intensity_matrix: np.ndarray, apex_idx: int,
                     n_avg_scans: int = 5, n_bg_scans: int = 3,
                     bg_offset_scans: int = 10) -> tuple:
    """Xcalibur Spectrum Enhancement: clean spectrum from a peak.

    1. Peak spectrum = mean of apex ± n_avg_scans//2 scans
    2. Left bg = mean of scans bg_offset_scans before peak
    3. Right bg = mean of scans bg_offset_scans after peak
    4. Background = (left_bg + right_bg) / 2
    5. Clean spectrum = peak_spectrum - background (clip to 0)

    Args:
        intensity_matrix: Full 2D matrix [n_scans, n_mz]
        apex_idx: Peak apex scan index
        n_avg_scans: Scans to average at apex (odd, 3-7)
        n_bg_scans: Baseline scans to average each side
        bg_offset_scans: Offset from apex for bg regions

    Returns:
        (enhanced_spectrum, peak_spectrum_raw, bg_spectrum)
        enhanced_spectrum normalized to base peak = 100
    """
    n_scans = intensity_matrix.shape[0]
    half_avg = n_avg_scans // 2

    # Peak top region
    p_start = max(0, apex_idx - half_avg)
    p_end = min(n_scans - 1, apex_idx + half_avg)
    peak_spec = intensity_matrix[p_start:p_end + 1, :].mean(axis=0)

    # Left baseline
    lb_center = max(0, apex_idx - bg_offset_scans)
    lb_start = max(0, lb_center - n_bg_scans // 2)
    lb_end = min(n_scans - 1, lb_start + n_bg_scans - 1)
    left_bg = intensity_matrix[lb_start:lb_end + 1, :].mean(axis=0)

    # Right baseline
    rb_center = min(n_scans - 1, apex_idx + bg_offset_scans)
    rb_start = max(0, rb_center - n_bg_scans // 2)
    rb_end = min(n_scans - 1, rb_start + n_bg_scans - 1)
    right_bg = intensity_matrix[rb_start:rb_end + 1, :].mean(axis=0)

    # Background = mean of both sides, scaled to 50% to avoid over-subtraction
    bg_spec = (left_bg + right_bg) / 2.0

    # Subtract background at 60% weight (conservative: preserves real compound ions)
    enhanced = peak_spec - bg_spec * 0.6
    enhanced = np.clip(enhanced, 0, None)

    # Normalize to base peak = 100
    if enhanced.max() > 0:
        enhanced_norm = enhanced / enhanced.max() * 100.0
    else:
        enhanced_norm = enhanced

    return enhanced_norm, peak_spec, bg_spec


def enhance_all_peaks(intensity_matrix: np.ndarray, peaks: list) -> list:
    """Apply spectrum enhancement to all detected peaks.

    Dynamically chooses n_avg_scans based on peak width.
    """
    for peak in peaks:
        peak_width = peak['end_idx'] - peak['start_idx']
        n_avg = min(max(3, peak_width // 3), 7)  # 3-7 scans
        n_bg = max(3, peak_width // 4)
        bg_off = max(10, peak_width)  # Further offset = cleaner background

        enhanced, raw, bg = enhance_spectrum(
            intensity_matrix, peak['apex_idx'],
            n_avg_scans=n_avg, n_bg_scans=n_bg,
            bg_offset_scans=bg_off
        )

        # Store active m/z ions (intensity > 1%)
        active = enhanced > 1.0
        peak['enhanced_mz'] = np.where(active)[0]  # bin indices
        peak['enhanced_int'] = enhanced[active]
        peak['enhanced_full'] = enhanced  # Full spectrum for later use
        peak['molecular_ion_candidate'] = (
            float(np.where(enhanced > 5.0)[0].max())
            if (enhanced > 5.0).any() else 0
        )

    return peaks
