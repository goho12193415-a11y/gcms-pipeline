"""
Step 8: Quantification
=======================
Three methods: external standard, internal standard, area normalization.
"""
import numpy as np


class QuantificationEngine:
    """GC-MS quantification engine."""

    @staticmethod
    def external_standard(sample_area: float,
                          calibration_curve: list) -> dict:
        """External standard: linear regression of standards.

        calibration_curve: [{level: mg/L, area: float}, ...]
        """
        if len(calibration_curve) < 2:
            raise ValueError("Need at least 2 calibration points")

        levels = np.array([p['level'] for p in calibration_curve])
        areas = np.array([p['area'] for p in calibration_curve])

        coeffs = np.polyfit(levels, areas, 1)
        k, b = coeffs
        y_pred = np.polyval(coeffs, levels)
        ss_res = np.sum((areas - y_pred) ** 2)
        ss_tot = np.sum((areas - np.mean(areas)) ** 2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0

        concentration = max(0.0, (sample_area - b) / k) if k != 0 else 0

        return {
            'concentration': concentration,
            'slope': k, 'intercept': b,
            'r_squared': r2,
            'method': 'external_standard'
        }

    @staticmethod
    def internal_standard(sample_area: float, is_area: float,
                          is_concentration: float,
                          response_factor: float) -> dict:
        """Internal standard: compensates injection volume & drift.

        C_analyte = (A_analyte / A_IS) × C_IS / RRF
        """
        if is_area == 0:
            return {
                'concentration': 0,
                'error': 'IS peak not found',
                'method': 'internal_standard'
            }

        area_ratio = sample_area / is_area
        concentration = (area_ratio * is_concentration) / response_factor

        return {
            'concentration': max(0, concentration),
            'area_ratio': area_ratio,
            'response_factor': response_factor,
            'method': 'internal_standard'
        }

    @staticmethod
    def area_normalization(target_area: float,
                           all_peak_areas: list) -> dict:
        """Area normalization: relative percentage of total area.

        Suitable for untargeted profiling without standards.
        """
        total = sum(all_peak_areas)
        if total == 0:
            return {'percentage': 0, 'method': 'normalization'}

        percentage = target_area / total * 100
        return {'percentage': percentage, 'method': 'normalization'}
