# GC-MS Auto-Identification Pipeline

Originally built for food and marine volatile organic compound (VOC) analysis — now expanded to 406 NIST-verified EI-MS spectra across 14 application domains: Maillard/thermal, meat, tea/coffee, spices, fruits, fermented beverages, dairy, lipid oxidation, nuts, smoke, vegetables, carotenoid degradation, and environmental contaminants.

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

A curated MSP-format spectral library of **300 compounds** organized by application domain. All spectra verified against NIST reference data for base peak accuracy.

### Quick Lookup by Analysis Type

| Application | Count | Key Markers |
|-------------|:-----:|-------------|
| **Marine & Algae** | 41 | DMS, halogenated phenols, iodoalkanes, 1-penten-3-ol, geosmin, (E,Z)-2,6-nonadienal |
| **Maillard & Thermal** | 65 | Pyrazines, furans, thiazoles, Strecker aldehydes, pyrroles, furanones |
| **Fruity & Floral** | 51 | Monoterpenes, ionones, damascenones, lactones, benzaldehyde, esters |
| **Lipid Oxidation** | 32 | Alkanals, 2-alkenals, 2,4-dienals, 1-octen-3-ol/one, alkylfurans |
| **Sesquiterpenes** | 30 | Caryophyllene, cadinenes, muuroles, elemenes, oxygenated sesquiterpenes |
| **Fermented & Spoilage** | 21 | Short-chain acids, ethyl esters, acetoin, 2,3-butanediol |
| **Environmental** | 17 | Phthalates, BHT, BTEX, naphthalenes |
| **Carotenoid Degradation** | 16 | Ionones, damascenones, safranal, megastigmatrienone |
| **Fresh & Green** | 15 | C6 aldehydes/alcohols, hexenals, hexenols |
| **Phenolic & Smoky** | 11 | Guaiacol, eugenol, vanillin, vinylguaiacol, cresols |

See [category_index.json](category_index.json) for the complete categorized compound list.

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
