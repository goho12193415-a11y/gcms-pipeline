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

# ---- Step 1: RAW input ----
# Read Thermo .RAW natively (Thermo RawFileReader via pythonnet/coreclr) instead
# of converting to mzML first. Same instrument data Xcalibur reads, one less
# dependency, faster. Verified byte-equivalent spectra to the mzML path.
NIST_NATIVE_RAW = True

# ---- Step 6a: NIST libraries ----
# Searching the replicate library (replib) alongside mainlib recovers compounds
# whose curated mainlib reference doesn't match this instrument but whose
# replicate spectrum does (verified +3 in ISOLATION on xt6.26: 2,3-octanedione,
# (Z)-7-hexadecenal, the diisobutyrate). BUT in a fixed top-5 review budget the
# +3 gains are offset by ~4 displacement losses (1-octen-3-one, PEG ethers,
# where mainlib already had the answer at rank 2-5): net 23->22. Tested with
# naive-merge, RRF, RMF-gated, and reserve-slot strategies — all net-negative
# at top-5. With the review window widened to top-8 (TOP_N_CANDIDATES) the
# replib recoveries land in rank 6-8 without displacing the mainlib hits, so it
# is a net gain — enabled together with TOP_N_CANDIDATES=8.
NIST_USE_REPLIB = True
NIST_REPLIB_CANDIDATES = [
    r"C:\Users\go ho\Desktop\MSSEARCH\replib",
    r"C:\Users\go ho\Desktop\GCMS_Software\谱库\replib",
]
# Only consult replib when mainlib is NOT already confident (top-1 RMF below
# this). Merging replib into confident mainlib peaks displaces correct hits
# (measured: 1-octen-3-one, PEG ethers). Gating captures the replib gains on
# weak peaks without the losses on strong ones.
NIST_REPLIB_FALLBACK_RMF = 850
NIST_REPLIB_RESERVE = 3     # slots reserved for replib-unique hits (keep=top_n-reserve mainlib ranks intact)
NIST_REPLIB_MIN_RMF = 750   # only inject a replib candidate if its RMF >= this (cuts review noise)

# ---- Step 6b: NIST query spectrum cleanup ----
# Air / background ions that dominate low-abundance peak spectra (base peak is
# water m/z 18, then N2/O2/Ar), burying the real characteristic ions and
# wrecking the NIST match. Stripped before normalization. Diagnosed as the root
# cause of the Bm-class failures (see 项目技术档案.md, 2026-06-29).
AIR_IONS = (17, 18, 28, 32, 40, 44)   # H2O, N2/CO, O2, Ar, CO2

# ---- Step 7: Library Matching ----
MIN_RMF = 700              # Minimum reverse match factor
MZ_WEIGHT_EXP = 3.0        # m/z weighting exponent (NIST standard = 3)
INT_WEIGHT_EXP = 0.6       # Intensity weighting exponent (NIST standard = 0.6)
RI_TOLERANCE = 30          # RI tolerance for hard filter
RI_WEIGHT = 0.20           # RI contribution (20% — RI is most reliable orthogonal evidence)

# ---- Step 7b: RI-primary injection for branched alkanes ----
# Branched-alkane homologs/isomers are spectrally near-identical and the exact
# compound is usually NOT in the NIST spectral top-50 (diagnosed: 0/11). Idea:
# when a peak's spectrum reads "alkane", inject the predicted-RI (CNN/Goodner)
# candidate(s) closest to the measured RI.
#
# VERDICT (2026-06-29): DISABLED. Tested rigorously, net zero on xt6.26 with
# false-positive risk. Two complementary, unfixable failure modes:
#   - where predicted RI is decisive (e.g. 2,6,10-Trimethylpentadecane |dRI|=4),
#     the apex spectrum does NOT read as an alkane (NIST returns alcohols), so
#     the alkane gate cannot fire;
#   - where the spectrum DOES read as alkane (2,6,10-Trimethyldodecane), the
#     correct homolog is buried among closer-RI competitors (isomer dRI < pred
#     error), so RI cannot pick it.
# Plus several A-class ground-truth labels are internally inconsistent (same
# name at impossible RIs). Code/flag kept for a future, more accurate RI model.
RI_INJECT_ALKANES = False  # RI-primary injection for branched alkanes (see above)
RI_INJECT_TOL = 25         # Max |predicted-RI - measured-RI| to inject (RI units)
RI_INJECT_MAX = 2          # Max candidates to inject per peak

# ---- Step 8: Quantification ----
QUANT_METHOD = "normalization"  # external, internal, or normalization

# ---- Step 9: Export ----
EXPORT_FORMAT = "xlsx"
# Number of NIST candidates retained per peak (search depth + Excel rows). The
# correct answer for weak/replib-recovered peaks can sit at rank 6-8, so the
# review window is widened from 3 to 8 to surface them.
TOP_N_CANDIDATES = 8

# ---- Step 10: QC ----
QC_SAMPLE_PATTERN = "QC"
RT_DRIFT_THRESHOLD = 0.05  # min
AREA_CV_THRESHOLD = 0.15   # 15% CV

# ---- RI Calibration (from zgwt alkane standard, carbon # verified by M+) ----
# Carbon numbers were confirmed directly from the molecular ions in zgwt
# (each peak has exactly one M+ = 14n+2): RT 4.84=C10(142), 7.20=C11(156),
# 9.86=C12(170)... The PREVIOUS table was off by one carbon (it labeled the
# 7.2 min peak C10), which made every measured RI ~100 too low. Corrected here.
ALKANE_RTS_DB_WAX = [
    (4.84, 1000), (7.20, 1100), (9.86, 1200), (12.56, 1300),
    (15.18, 1400), (17.69, 1500), (20.08, 1600), (22.36, 1700),
    (24.54, 1800), (26.79, 1900), (29.03, 2000), (31.12, 2100),
    (33.08, 2200), (34.99, 2300), (37.32, 2400), (39.43, 2500),
    (41.09, 2600), (42.50, 2700),
]
RI_COLUMN = "DB-WAX"

# ---- RI QC: this-column <-> literature DB-WAX transform ----
# Literature DB-WAX RI does NOT transfer 1:1 to this column (different RI scale).
# Fitted on 15 distinctive high-RMF compounds (R^2=0.994):
#     literature_RI = RI_QC_SLOPE * measured_RI + RI_QC_INTERCEPT
# A correct ID lies on this line; |residual| > RI_QC_TOL flags a suspect ID.
RI_QC_SLOPE = 0.635
RI_QC_INTERCEPT = 490
RI_QC_TOL = 30
