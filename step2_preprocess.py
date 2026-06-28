"""
Step 2: Signal Preprocessing (v2 - SNIP + S-G + Column Bleed)
===============================================================
Per GCMS_Complete_Workflow_v2.md Layer 0, 2:
  - Column bleed spectrum estimation & subtraction
  - Savitzky-Golay smoothing (vectorized)
  - SNIP baseline correction (vectorized, ~50× faster than AsLS)
"""
import numpy as np
from scipy.signal import savgol_filter


# ---- Column Bleed (Layer 0.4) ----

COLUMN_BLEED_MZ = {73, 147, 207, 221, 281, 355, 429}


def estimate_column_bleed(intensity_matrix: np.ndarray,
                           late_rt_fraction: float = 0.1) -> np.ndarray:
    """Estimate column bleed spectrum from late-eluting scans.

    Takes mean spectrum of the last 10% of scans (highest temperature,
    most column bleed). This is the Xcalibur auto-background logic.
    """
    n_scans = intensity_matrix.shape[0]
    late_start = int(n_scans * (1 - late_rt_fraction))
    return intensity_matrix[late_start:, :].mean(axis=0)


def subtract_column_bleed(intensity_matrix: np.ndarray,
                           bleed_spectrum: np.ndarray,
                           scale: float = 0.9) -> np.ndarray:
    """Subtract scaled column bleed spectrum from all scans."""
    corrected = intensity_matrix - bleed_spectrum[np.newaxis, :] * scale
    return np.clip(corrected, 0, None)


# ---- Savitzky-Golay Smoothing (Layer 2.1) ----

def smooth_signal(signal: np.ndarray, window_length: int = 5,
                  polyorder: int = 2) -> np.ndarray:
    """S-G filter: retains peak shape, suppresses high-freq noise."""
    if len(signal) < window_length:
        return signal
    return np.clip(savgol_filter(signal, window_length, polyorder), 0, None)


# ---- SNIP Baseline Correction (Layer 2.2) ----

def snip_baseline(signal: np.ndarray, max_half_window: int = 100) -> np.ndarray:
    """SNIP: Statistics-sensitive Non-linear Iterative Peak-clipping.

    Vectorized implementation. For each iteration p:
      neighbor_mean = (y[i-p] + y[i+p]) / 2
      if y[i] > neighbor_mean: y[i] = neighbor_mean

    This progressively clips peaks, leaving the true baseline.
    ~50× faster than AsLS sparse solver.
    """
    y = np.copy(signal).astype(float)
    n = len(y)

    for p in range(1, min(max_half_window + 1, n // 2)):
        y_left = np.empty(n)
        y_right = np.empty(n)
        y_left[:p] = y[:p]
        y_right[-p:] = y[-p:]
        y_left[p:] = y[:-p]
        y_right[:-p] = y[p:]
        neighbor_mean = (y_left + y_right) / 2.0
        y = np.minimum(y, neighbor_mean)

    return y


def correct_baseline(signal: np.ndarray,
                      max_half_window: int = 100) -> tuple:
    """Full SNIP baseline correction pipeline."""
    baseline = snip_baseline(signal, max_half_window)
    corrected = np.clip(signal - baseline, 0, None)
    return corrected, baseline


def preprocess_signal(signal: np.ndarray, window_length: int = 5,
                      polyorder: int = 2,
                      max_half_window: int = 100) -> dict:
    """Run full Step 2 preprocessing (S-G smooth → SNIP baseline)."""
    smoothed = smooth_signal(signal, window_length, polyorder)
    corrected, baseline = correct_baseline(smoothed, max_half_window)
    return {
        'smoothed': smoothed,
        'baseline': baseline,
        'corrected': corrected
    }
