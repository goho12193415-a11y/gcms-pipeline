#!/usr/bin/env python3
"""
Convert NIST SDF spectral library to MSP format for use with gcms_pipeline.py.

Usage:
    python convert_nist_sdf.py mainlib.sdf nist_mainlib.msp

The NIST 2014 (or later) EI MS library in SDF format can be found in the
Supplemental.zip included with NIST MS Search installations.

This extracts all organic compounds (MW < 600, containing carbon) and
writes them in MSP format compatible with gcms_pipeline.py.

WARNING: The NIST library is copyrighted. Do not redistribute the
converted MSP file. Use only for personal/academic research.
"""

import re, sys, os

def convert_sdf_to_msp(sdf_path, msp_path, max_mw=600):
    print(f"Reading: {sdf_path}")
    with open(sdf_path, 'r', encoding='latin-1', errors='replace') as f:
        text = f.read()

    entries = text.split('$$$$\r\n')
    print(f"Entries found: {len(entries)}")

    converted = 0
    skipped = 0

    with open(msp_path, 'w', encoding='utf-8') as fout:
        for i, entry in enumerate(entries):
            if not entry.strip():
                continue
            if i % 50000 == 0 and i > 0:
                print(f"  Processing {i:,}/{len(entries)}... ({converted:,} converted)")

            name_m = re.search(r'>  <NAME>\r?\n(.+)', entry)
            cas_m = re.search(r'>  <CASNO>\r?\n(.+)', entry)
            formula_m = re.search(r'>  <FORMULA>\r?\n(.+)', entry)
            mw_m = re.search(r'>  <MW>\r?\n(.+)', entry)
            peaks_m = re.search(
                r'>  <MASS SPECTRAL PEAKS>\r?\n(.*?)(?:\r?\n\r?\n|\r?\n$)', entry, re.DOTALL
            )

            if not name_m or not peaks_m:
                skipped += 1; continue

            name = name_m.group(1).strip()
            formula = formula_m.group(1).strip() if formula_m else ''
            cas = cas_m.group(1).strip() if cas_m else ''
            try:
                mw = float(mw_m.group(1)) if mw_m else 999
            except:
                mw = 999

            # Filter: MW limit, require carbon (organic)
            if mw > max_mw or 'C' not in formula:
                skipped += 1; continue

            # Parse peaks
            peak_lines = peaks_m.group(1).strip().split('\n')
            peaks = []
            for pl in peak_lines:
                parts = pl.strip().split()
                if len(parts) >= 2:
                    try:
                        mz = int(parts[0]); intens = int(parts[1])
                        if mz > 0 and intens > 0:
                            peaks.append((mz, intens))
                    except:
                        pass

            if len(peaks) < 3:
                skipped += 1; continue

            # Normalize and keep top 50 peaks
            max_int = max(i for _, i in peaks)
            peaks_norm = [(m, int(i / max_int * 999)) for m, i in peaks
                         if int(i / max_int * 999) > 0]
            top_peaks = sorted(peaks_norm, key=lambda x: -x[1])[:50]

            fout.write(f"Name: {name}\n")
            fout.write(f"Formula: {formula}\n")
            fout.write(f"CASNO: {cas}\n")
            fout.write(f"Num Peaks: {len(top_peaks)}\n")
            fout.write('; '.join(f"{m} {i}" for m, i in top_peaks) + '\n')
            converted += 1

    sz = os.path.getsize(msp_path) / 1024 / 1024
    print(f"\nDone! Converted: {converted:,}, Skipped: {skipped:,}")
    print(f"Output: {msp_path} ({sz:.0f} MB)")


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(__doc__)
        print("Usage: python convert_nist_sdf.py <mainlib.sdf> [output.msp]")
        sys.exit(1)

    sdf_in = sys.argv[1]
    msp_out = sys.argv[2] if len(sys.argv) > 2 else 'nist_mainlib.msp'
    convert_sdf_to_msp(sdf_in, msp_out)
