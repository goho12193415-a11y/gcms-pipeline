"""
GC-MS Pipeline v2 — Configuration
==================================
All tunable parameters with defaults from the workflow document.
"""
from pathlib import Path

# ---- Paths ----
PROJECT_DIR = Path(__file__).parent
LIBRARY_PATH = PROJECT_DIR / "library" / "nist_gcms.json"  # Full NIST (239K)
MZML_DIR = Path(r"C:\Users\go ho\Desktop\gcms_mzml")
OUTPUT_DIR = PROJECT_DIR / "output"

# ---- Step 2: Signal Preprocessing ----
SMOOTH_WINDOW = 5          # Savitzky-Golay window (must be odd)
SMOOTH_POLYORDER = 2       # Polynomial order for S-G filter
BASELINE_LAMBDA = 1e5      # AsLS smoothness (1e3~1e7, larger = flatter)
BASELINE_P = 0.001         # AsLS asymmetry (0.001~0.01)

# ---- Step 3: Peak Detection ----
MIN_SN = 10.0              # Minimum signal-to-noise ratio (was 5.0)
MIN_PEAK_WIDTH_SCANS = 3   # Minimum peak width in scan points
SMOOTH_SIGMA = 1.0         # Gaussian smoothing sigma before peak detection

# ---- Step 4: Deconvolution ----
MZ_MIN = 35                # Minimum m/z for intensity matrix
MZ_MAX = 500               # Maximum m/z for intensity matrix
MZ_STEP = 1.0              # m/z bin width
CORRELATION_THRESHOLD = 0.95  # Ion EIC correlation for clustering (was 0.80)
DECONV_ENABLED = False        # Deconvolution off by default (--deconv to enable)
DECONV_VALLEY_RATIO = 0.8     # Only deconvolve if valley/peak ratio > 0.8 (co-elution)
MIN_IONS = 3               # Minimum ions per deconvoluted component

# ---- Step 5: Integration ----
USE_PEAK_BASELINE = True   # Local linear baseline per peak
INTEGRATION_METHOD = "trapezoid"  # trapezoid or simpson

# ---- Step 6: RT Alignment ----
RT_ALIGN_TOLERANCE = 0.2   # Initial matching window (min)
RT_ALIGN_POLY_DEGREE = 1   # 1=linear, 2=quadratic

# ---- Step 7: Library Matching ----
MIN_RMF = 700              # Minimum reverse match factor
MZ_WEIGHT_EXP = 3.0        # m/z weighting exponent (NIST standard = 3)
INT_WEIGHT_EXP = 0.6       # Intensity weighting exponent (NIST standard = 0.6)
RI_TOLERANCE = 30          # RI tolerance for hard filter
RI_WEIGHT = 0.20           # RI contribution (20% — RI is most reliable orthogonal evidence)

# ---- Step 8: Quantification ----
QUANT_METHOD = "normalization"  # external, internal, or normalization

# ---- Step 9: Export ----
EXPORT_FORMAT = "xlsx"

# ---- Step 10: QC ----
QC_SAMPLE_PATTERN = "QC"
RT_DRIFT_THRESHOLD = 0.05  # min
AREA_CV_THRESHOLD = 0.15   # 15% CV

# ---- RI Calibration (from zgwt alkane standard) ----
# DB-WAX, 40→230°C, calibrated via Benzaldehyde (Δ=19) & Nonanal (Δ=74)
ALKANE_RTS_DB_WAX = [
    (7.2, 1000), (9.9, 1100), (12.6, 1200), (15.2, 1300),
    (17.7, 1400), (20.1, 1500), (22.4, 1600), (24.5, 1700),
    (26.8, 1800), (29.0, 1900), (31.1, 2000), (33.1, 2100),
    (35.0, 2200), (37.3, 2300), (39.4, 2400), (41.1, 2500),
    (42.5, 2600), (44.0, 2700),
]
RI_COLUMN = "DB-WAX"
