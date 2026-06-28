"""
Step 5: Peak Area Integration
==============================
Trapezoidal rule integration with local linear baseline per peak.
"""
import numpy as np


def integrate_peak(
    signal: np.ndarray,
    rt: np.ndarray,
    peak: dict,
    use_peak_baseline: bool = True
) -> dict:
    """Integrate a single chromatographic peak.

    Connects start→end points as local linear baseline,
    subtracts it, then applies trapezoidal integration.

    Args:
        signal: Baseline-corrected TIC signal
        rt: Retention time array
        peak: Peak dict from Step 3
        use_peak_baseline: If True, draw local baseline between
                          peak start and end points

    Returns:
        Peak dict with added 'area' and 'height' keys
    """
    start = peak['start_idx']
    end = peak['end_idx'] + 1

    peak_rt = rt[start:end]
    peak_signal = signal[start:end]

    if len(peak_rt) < 2:
        return {**peak, 'area': 0.0,
                'height': float(signal[peak['apex_idx']])}

    if use_peak_baseline:
        y_start = float(peak_signal[0])
        y_end = float(peak_signal[-1])
        baseline_line = np.linspace(y_start, y_end, len(peak_signal))
        net_signal = np.clip(peak_signal - baseline_line, 0, None)
    else:
        net_signal = peak_signal

    area = float(np.trapz(net_signal, peak_rt))
    height = float(signal[peak['apex_idx']])

    return {**peak, 'area': area, 'height': height}


def integrate_all_peaks(signal: np.ndarray, rt: np.ndarray,
                         peaks: list, **kwargs) -> list:
    """Integrate all detected peaks."""
    return [integrate_peak(signal, rt, p, **kwargs) for p in peaks]
