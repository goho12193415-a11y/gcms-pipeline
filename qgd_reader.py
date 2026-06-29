"""
qgd_reader.py — native Shimadzu GCMSsolution .qgd reader.
=========================================================
.qgd is an OLE2 compound file. The MS scan data lives in:
  GCMS Raw Data/Retention Time  : int32 ms per scan
  GCMS Raw Data/Spectrum Index  : int32 byte offset per scan into MS Raw Data
  GCMS Raw Data/MS Raw Data     : per-scan [header][ (uint16 m/z*20, uintN intensity) ... ]
  GCMS Raw Data/TIC Data        : int32 TIC per scan (every other word) — used as
                                  a per-scan CHECKSUM to pick the right decoding.

Intensity is stored adaptively: 2 bytes (header 24) for low-intensity scans,
3 bytes (header 32) for high-intensity scans. For each scan we try both and
accept the variant whose intensity sum equals the stored TIC exactly — so a
mis-decode can never pass silently. Returns the same structure as
step1_parse.load_mzml_to_matrix.
"""
from pathlib import Path
import numpy as np
import olefile

_VARIANTS = [(24, 4), (32, 5)]   # (header_bytes, point_bytes); point = u16 mz*20 + uintN

# Shimadzu raw scans carry a flat low-level noise floor (a roughly CONSTANT
# ~500-count value at every m/z channel), unlike Thermo mzML which is
# pre-thresholded. Left in, those hundreds of spurious ions wreck the NIST
# reverse-match (it matches big molecules with ions everywhere). The floor is
# constant in absolute terms, so a per-scan median captures it well across all
# abundances; keep points above NOISE_MULT * median. Recovers weak-peak IDs
# (Octanal/1-Hexanol) that a relative-to-base threshold missed.
_NOISE_MULT = 3.0


def _parse(blk, H, w):
    body = blk[H:]
    m = len(body) // w
    if m == 0:
        return np.empty(0, np.float32), np.empty(0, np.float32)
    a = np.frombuffer(body[:m * w], np.uint8).reshape(-1, w)
    mz = (a[:, 0].astype(np.uint32) | a[:, 1].astype(np.uint32) << 8) / 20.0
    it = a[:, 2].astype(np.uint32) | a[:, 3].astype(np.uint32) << 8
    if w == 5:
        it |= a[:, 4].astype(np.uint32) << 16
    return mz.astype(np.float32), it.astype(np.float32)


def load_qgd_to_matrix(path: str) -> dict:
    ole = olefile.OleFileIO(str(path))
    rt_ms = np.frombuffer(ole.openstream('GCMS Raw Data/Retention Time').read(), '<i4')
    idx = np.frombuffer(ole.openstream('GCMS Raw Data/Spectrum Index').read(), '<i4')
    raw = ole.openstream('GCMS Raw Data/MS Raw Data').read()
    tic = np.frombuffer(ole.openstream('GCMS Raw Data/TIC Data').read(), '<i4')[0::2]
    ole.close()

    N = len(idx)
    scan_list, rt_list, tic_list = [], [], []
    n_unverified = 0
    for i in range(N):
        s = idx[i]
        e = idx[i + 1] if i + 1 < N else len(raw)
        blk = raw[s:e]
        chosen = None
        for (H, w) in _VARIANTS:
            mz, it = _parse(blk, H, w)
            if abs(float(it.sum()) - float(tic[i])) < 2:
                chosen = (mz, it)
                break
        if chosen is None:                       # TIC checksum failed — flag, best-effort
            n_unverified += 1
            chosen = _parse(blk, 24, 4)
        mz, it = chosen
        thr = _NOISE_MULT * np.median(it) if it.size else 0.0   # strip noise floor
        k = it > thr
        mz, it = mz[k], it[k]
        rt = float(rt_ms[i]) / 60000.0           # ms -> min
        rt_list.append(rt)
        tic_list.append(float(it.sum()))
        scan_list.append({'rt': rt, 'mz': mz, 'intensity': it})

    if n_unverified:
        print(f"  [qgd] WARNING: {n_unverified}/{N} scans failed TIC checksum")
    mz_all = np.concatenate([s['mz'] for s in scan_list if len(s['mz']) > 0])
    return {
        'rt': np.array(rt_list, np.float64),
        'tic': np.array(tic_list, np.float64),
        'scan_list': scan_list,
        'mz_range': (float(mz_all.min()), float(mz_all.max())),
        'sample_name': Path(path).stem,
        'n_scans': N,
    }


if __name__ == "__main__":
    import sys
    d = load_qgd_to_matrix(sys.argv[1])
    print(f"{d['n_scans']} scans, RT {d['rt'][0]:.2f}-{d['rt'][-1]:.2f} min, "
          f"m/z {d['mz_range'][0]:.0f}-{d['mz_range'][1]:.0f}")
