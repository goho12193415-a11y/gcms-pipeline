# GC-MS Auto-Identification Pipeline

Automated peak detection, spectral matching, and compound identification for GC-MS data. Built for food and marine volatile organic compound (VOC) analysis.

## Features

- **Raw Data Parsing** — Reads Thermo Xcalibur GC-MS text exports (scan-level m/z + intensity)
- **Peak Detection & Integration** — Savitzky-Golay smoothing, baseline correction, Simpson integration
- **Spectral Library Matching** — NIST-style sqrt-weighted cosine similarity with forward/reverse search
- **Background Filtering** — Air/water, column bleed, phthalate, and contaminant removal
- **RT Sanity Check** — Retention time validation for physically impossible matches
- **Excel Reporting** — Formatted reports with compound name, CAS, SI score, peak area, confidence level
- **Batch Processing** — Process multiple files in one command

## Installation

```bash
pip install numpy scipy openpyxl
```

## Quick Start

```bash
# Single file
python gcms_pipeline.py sample.txt output.xlsx

# Batch processing
python gcms_pipeline.py --batch "data/*.txt" ./results/

# With custom SI threshold (default: 850)
python gcms_pipeline.py sample.txt output.xlsx --threshold 800
```

## Input Format

Thermo Xcalibur GC-MS text export with scan headers and data packets:

```
ScanHeader # 1
start_time = 3.000832
integ_intens = 17755716.754425
...
Packet # 0, intensity = 821075.187500, mass/position = 17.118618
Packet # 1, intensity = 3683056.250000, mass/position = 18.065636
...
```

To export from Xcalibur: `File → Export → Text (.txt)` with "All scans" selected.

## Spectral Library

`food_volatiles.msp` contains **200 compounds** with NIST-verified reference EI-MS spectra:

| Class | Count | Examples |
|-------|:-----:|----------|
| Aldehydes | 25 | Hexanal, (E)-2-Hexenal, Nonanal, (E,Z)-2,6-Nonadienal |
| Ketones | 20 | 2-Heptanone, beta-Ionone, 6-Methyl-5-hepten-2-one |
| Alcohols | 18 | 1-Octen-3-ol, 1-Penten-3-ol, Geosmin |
| Terpenes | 55 | Limonene, beta-Caryophyllene, Germacrene D, Phytol |
| Sulfur Compounds | 8 | Dimethyl sulfide, Dimethyl disulfide, Benzothiazole |
| Halogenated | 12 | Iodoform, 1-Iodoheptane, 2,6-Dibromophenol |
| Pyrazines | 8 | 2,5-Dimethylpyrazine, 2,3,5-Trimethylpyrazine |
| Furans | 10 | 2-Pentylfuran, Furfural, 2-Acetylfuran |
| Phenols | 12 | Guaiacol, Eugenol, 4-Vinylguaiacol |
| Acids | 12 | Acetic acid, Hexanoic acid, Nonanoic acid |
| Esters/Lactones | 10 | gamma-Butyrolactone, gamma-Nonalactone |
| N-Compounds | 6 | Indole, 2-Acetylpyrrole, 3-Methylindole |
| Others | 9 | BHT, 2-MIB, Neophytadiene, Squalene |

### Expanding the Library

To add new compounds, append to `food_volatiles.msp` in MSP format:

```
Name: Compound Name
Formula: C6H12O
CASNO: 66-25-1
Num Peaks: 18
44 999; 56 450; 41 400; 43 350; 72 300
```

Peak format: `m/z relative_intensity` pairs separated by `;`.

## Output

Excel report with columns:
- **RT (min)** — Retention time
- **Compound** — Best library match
- **CAS No.** — CAS registry number
- **SI** — Spectral match factor (0-999, NIST-style)
- **Confidence** — High (SI≥880) / Medium (SI≥850) / Low (SI<850)
- **Peak Area** — Integrated peak area
- **Peak Height** — TIC at peak apex

## Limitations

- Library quality determines identification accuracy. Compounds not in the library cannot be identified.
- Early retention time (RT < 6 min) and late (RT > 40 min) regions are prone to false positives from air/water and column artifacts.
- RT-based filtering relies on DB-5 MS column with ~50 min temperature program. Adjust `RT_RANGES` for different conditions.
- For publication-quality identification, verify key compounds against NIST library in Xcalibur.

## Citation

If you use this pipeline in your research, please cite:
```
GC-MS Auto-Identification Pipeline. https://github.com/goho12193415-a11y/gcms-pipeline
```

Author contact: goho12193415@gmail.com

## License

MIT License. The bundled `food_volatiles.msp` spectral library is freely redistributable.
