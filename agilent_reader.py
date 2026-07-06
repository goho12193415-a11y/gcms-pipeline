#!/usr/bin/env python3
"""
agilent_reader.py — native reader for Agilent MassHunter GC-MS/GC-Q-TOF .D data.
================================================================================
Reads an Agilent ``.D`` acquisition folder (MassHunter format) and returns the
same dict structure as qgd_reader / raw_reader, so the pipeline can process
Agilent data alongside Thermo .RAW and Shimadzu .qgd.

Format (reverse-engineered, validated against the values MassHunter stores):
  AcqData/MSScan.bin   — one ScanRecord per scan (RT, TIC, BasePeakMZ, and a
                         SpectrumParamValues block list). Parsed via rainbow's
                         scan-record reader (the schema is in MSScan.xsd).
  AcqData/MSPeak.bin    — per scan, the centroid block (SpectrumFormatID == 2):
                         PointCount doubles of a raw time-of-flight axis X,
                         followed by PointCount float32 intensities.
  AcqData/DefaultMassCal.xml — the traditional TOF calibration:
                         m/z = (coeff * (X - base))**2      (Step 1, 2 coeffs)

Verified on SPME GC-Q-TOF data: reconstructed TIC and base-peak m/z match the
values MassHunter itself recorded in each ScanRecord exactly. m/z is binned to
nominal mass (unit) for NIST EI search, matching the other readers.
"""
import os
from pathlib import Path
import numpy as np
from xml.etree import ElementTree as ET

from rainbow.agilent import masshunter as _mh  # scan-index parsing (schema-driven)


def _find_acqdata(path):
    """Accept either the .D folder or the AcqData folder itself."""
    p = Path(path)
    if (p / "AcqData" / "MSScan.bin").exists():
        return str(p / "AcqData")
    if (p / "MSScan.bin").exists():
        return str(p)
    raise FileNotFoundError(f"No AcqData/MSScan.bin under {path}")


def _read_traditional_cal(acqdata):
    """(coeff, base) of the traditional TOF calibration from DefaultMassCal.xml.
    Falls back to identity-ish if absent (then X is treated as m/z)."""
    xml = os.path.join(acqdata, "DefaultMassCal.xml")
    if not os.path.exists(xml):
        return None
    root = ET.parse(xml).getroot()
    # first DefaultCalibration -> Step with Traditional formula -> 2 Values
    for cal in root.iter("DefaultCalibration"):
        for step in cal.iter("Step"):
            if step.findtext("CalibrationFormula") == "Traditional":
                vals = [float(v.text) for v in step.iter("Value")]
                if len(vals) >= 2:
                    return vals[0], vals[1]      # coeff, base
    return None


def load_agilent_to_matrix(path: str) -> dict:
    acqdata = _find_acqdata(path)
    cal = _read_traditional_cal(acqdata)

    ctypes = _mh.parse_scan_xsd(os.path.join(acqdata, "MSScan.xsd"))
    records = _mh.read_scan_records(
        os.path.join(acqdata, "MSScan.bin"), ctypes, _mh.count_scans(acqdata))
    with open(os.path.join(acqdata, "MSPeak.bin"), "rb") as f:
        peak_bytes = f.read()
    peak_size = len(peak_bytes)

    def to_mz(X):
        if cal is None:
            return np.asarray(X, np.float64)
        coeff, base = cal
        return (coeff * (np.asarray(X, np.float64) - base)) ** 2

    rt_list, tic_list, scan_list = [], [], []
    n_bad = 0
    for rec in records:
        rt = float(rec["ScanTime"])                 # minutes
        block = None
        for b in rec.get("SpectrumParamsBlocks", []):
            if b.get("SpectrumFormatID") == 2:      # centroid block
                block = b
                break
        if block is None:
            n_bad += 1
            rt_list.append(rt); tic_list.append(0.0)
            scan_list.append({'rt': rt, 'mz': np.array([]), 'intensity': np.array([])})
            continue

        npk = int(block["PointCount"])
        off = int(block["SpectrumOffset"])
        bc = int(block["ByteCount"])
        if npk <= 0 or bc <= 0 or off + bc > peak_size or bc // npk != 12:
            n_bad += 1
            rt_list.append(rt); tic_list.append(0.0)
            scan_list.append({'rt': rt, 'mz': np.array([]), 'intensity': np.array([])})
            continue

        raw = peak_bytes[off:off + bc]
        X = np.frombuffer(raw[:npk * 8], "<f8")
        inten = np.frombuffer(raw[npk * 8:npk * 12], "<f4").astype(np.float64)
        mz = to_mz(X)

        # bin to nominal mass (unit resolution) for NIST EI search
        nom = np.rint(mz).astype(np.int64)
        keep = (nom >= 1) & (inten > 0)
        nom, iv = nom[keep], inten[keep]
        if nom.size:
            order = np.argsort(nom, kind="stable")
            nom, iv = nom[order], iv[order]
            uniq, start = np.unique(nom, return_index=True)
            summed = np.add.reduceat(iv, start)
            mz_arr = uniq.astype(np.float64)
            it_arr = summed
        else:
            mz_arr = np.array([]); it_arr = np.array([])

        rt_list.append(rt)
        tic_list.append(float(it_arr.sum()))
        scan_list.append({'rt': rt, 'mz': mz_arr, 'intensity': it_arr})

    if n_bad:
        print(f"  [agilent] {n_bad}/{len(records)} scans had no centroid block")
    nonempty = [s['mz'] for s in scan_list if len(s['mz']) > 0]
    mz_all = np.concatenate(nonempty) if nonempty else np.array([0.0])
    return {
        'rt': np.array(rt_list, np.float64),
        'tic': np.array(tic_list, np.float64),
        'scan_list': scan_list,
        'mz_range': (float(mz_all.min()), float(mz_all.max())),
        'sample_name': Path(path).name.replace('.D', '').replace('.d', ''),
        'n_scans': len(records),
    }


if __name__ == "__main__":
    import sys
    d = load_agilent_to_matrix(sys.argv[1])
    print(f"{d['sample_name']}: {d['n_scans']} scans, "
          f"RT {d['rt'][0]:.2f}-{d['rt'][-1]:.2f} min, "
          f"m/z {d['mz_range'][0]:.0f}-{d['mz_range'][1]:.0f}")
