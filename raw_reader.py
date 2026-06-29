"""
raw_reader.py — native Thermo .RAW reader (no mzML conversion, no Xcalibur).
============================================================================
Reads Thermo .RAW directly via Thermo's own RawFileReader .NET library
(bundled with ThermoRawFileParser) through pythonnet on the CoreCLR runtime.
This is the same instrument-native data Xcalibur/Qual Browser reads, with no
RAW->mzML conversion step. Returns the exact structure of
step1_parse.load_mzml_to_matrix so it is a drop-in replacement.

NOTE: pythonnet must load the CoreCLR runtime BEFORE `import clr`, because the
RawFileReader assemblies (v8) use C# default-interface-methods that the .NET
Framework runtime cannot load.
"""
from pathlib import Path
import numpy as np

_DLL_NAME = "ThermoFisher.CommonCore.RawFileReader.dll"
# Searched in order; first dir that contains the DLL wins. The bundled
# 'thermo_lib' next to this file makes the distribution self-contained; the
# C:\Tools path is the dev/ThermoRawFileParser fallback.
_DLL_DIR_CANDIDATES = [
    Path(__file__).resolve().parent / "thermo_lib",
    Path(r"C:\Tools\ThermoRawFileParser"),
]
_READY = False
_RawFileReaderAdapter = None
_Device = None


def _find_dll_dir():
    cands = list(_DLL_DIR_CANDIDATES)
    try:
        from config import RAW_READER_DLL_DIRS
        cands = [Path(d) for d in RAW_READER_DLL_DIRS] + cands
    except Exception:
        pass
    for d in cands:
        if (d / _DLL_NAME).exists():
            return d
    raise RuntimeError(
        "Thermo RawFileReader DLLs not found (looked in: "
        + "; ".join(str(d) for d in cands) + "). Put "
        "ThermoFisher.CommonCore.Data.dll + RawFileReader.dll in a 'thermo_lib' "
        "folder next to raw_reader.py, or set NIST_NATIVE_RAW=False to use mzML.")


def _ensure_clr():
    global _READY, _RawFileReaderAdapter, _Device
    if _READY:
        return
    d = _find_dll_dir()
    from pythonnet import load
    load("coreclr")
    import clr
    clr.AddReference(str(d / "ThermoFisher.CommonCore.Data.dll"))
    clr.AddReference(str(d / _DLL_NAME))
    from ThermoFisher.CommonCore.RawFileReader import RawFileReaderAdapter
    from ThermoFisher.CommonCore.Data.Business import Device
    _RawFileReaderAdapter = RawFileReaderAdapter
    _Device = Device
    _READY = True


def load_raw_to_matrix(raw_path: str) -> dict:
    """Read a Thermo .RAW natively. Same return shape as load_mzml_to_matrix."""
    _ensure_clr()
    raw = _RawFileReaderAdapter.FileFactory(str(raw_path))
    if raw.IsError:
        raise RuntimeError(f"RawFileReader cannot open {raw_path}")
    raw.SelectInstrument(_Device.MS, 1)
    hdr = raw.RunHeaderEx
    first, last = hdr.FirstSpectrum, hdr.LastSpectrum

    rt_list, tic_list, scan_list = [], [], []
    for sc in range(first, last + 1):
        rt = float(raw.RetentionTimeFromScanNumber(sc))
        stats = raw.GetScanStatsForScanNumber(sc)
        seg = raw.GetSegmentedScanFromScanNumber(sc, stats)
        mz = np.fromiter(seg.Positions, dtype=np.float32)
        inten = np.fromiter(seg.Intensities, dtype=np.float32)
        rt_list.append(rt)
        tic_list.append(float(inten.sum()))
        scan_list.append({'rt': rt, 'mz': mz, 'intensity': inten})
    raw.Dispose()

    mz_all = np.concatenate([s['mz'] for s in scan_list if len(s['mz']) > 0])
    return {
        'rt': np.array(rt_list, dtype=np.float64),
        'tic': np.array(tic_list, dtype=np.float64),
        'scan_list': scan_list,
        'mz_range': (float(mz_all.min()), float(mz_all.max())),
        'sample_name': Path(raw_path).stem,
        'n_scans': len(scan_list),
    }


if __name__ == "__main__":
    import sys, time
    t0 = time.time()
    d = load_raw_to_matrix(sys.argv[1])
    print(f"loaded {d['n_scans']} scans, RT {d['rt'][0]:.1f}-{d['rt'][-1]:.1f} min, "
          f"m/z {d['mz_range'][0]:.1f}-{d['mz_range'][1]:.1f} in {time.time()-t0:.1f}s")
