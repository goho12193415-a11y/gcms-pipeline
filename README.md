# GC-MS Auto-Identification Pipeline v2.1

Automated GC-MS data processing: `.RAW` → color-coded Excel report.

**Engine**: NIST MS Search (official DLL) — Chromeleon-grade spectral matching.
**Pipeline**: SNIP baseline → ICIS peaks → Spectrum enhancement → NIST search → Excel.

## Quick Start
```bash
pip install numpy scipy pandas openpyxl pymzml pyms-nist-search
python pipeline.py --mzml-files sample.mzML --nist
```

## What It Does
- RAW → mzML auto-conversion
- SNIP baseline correction (vectorized)
- ICIS peak detection (385 peaks from 13K scans)
- Spectrum enhancement (multi-scan average − background subtraction)
- NIST MS Search engine (official DLL, raw m/z precision)
- Kovats RI calculation (C10-C27 alkane calibration)
- Dual-column RI annotation (DB-WAX + DB-5)
- Review-optimized Excel (GREEN/YELLOW/RED/GRAY color coding)
- ~15 sec per sample

## Limitations
- Branched alkanes: EI-MS physical limit (ACS Anal. Chem. 2016)
- Automated ID ceiling: ~27% for complex food volatiles
- Pipeline is a productivity tool — expert review required for final ID

## License
MIT
