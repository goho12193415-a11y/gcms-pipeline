"""
Step 1: RAW → mzML conversion + data loading
=============================================
Uses ThermoRawFileParser for conversion, pymzml for loading.
"""
import subprocess
import numpy as np
from pathlib import Path
import pymzml


THERMO_PARSER_DLL = r"C:\Tools\ThermoRawFileParser\ThermoRawFileParser.dll"


def convert_raw_to_mzml(raw_path: str, output_dir: str) -> str:
    """Convert Thermo .RAW to .mzML via ThermoRawFileParser.

    Args:
        raw_path: Path to .RAW file
        output_dir: Directory for output .mzML

    Returns:
        Path to the generated .mzML file
    """
    raw_path = Path(raw_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "dotnet", THERMO_PARSER_DLL,
        f"-i={raw_path}",
        f"-o={output_dir}",
        "-f=1",          # 1 = mzML format
        "-p"             # Disable Thermo native peak picking (keep profile data)
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"RAW→mzML conversion failed:\n{result.stderr}")

    mzml_path = output_dir / (raw_path.stem + ".mzML")
    if not mzml_path.exists():
        raise FileNotFoundError(f"mzML not generated at {mzml_path}")

    return str(mzml_path)


def load_mzml_to_matrix(mzml_path: str) -> dict:
    """Read mzML file and build in-memory data matrices.

    Returns:
        data = {
            'rt': np.array,           # retention time series (min)
            'tic': np.array,          # TIC intensity series
            'scan_list': list[dict],  # each scan: {rt, mz: np.array, intensity: np.array}
            'mz_range': (min, max),   # m/z range
            'sample_name': str,
            'n_scans': int
        }
    """
    run = pymzml.run.Reader(mzml_path, MS_precision={1: 5e-3},
                           build_index_from_scratch=True)

    rt_list, tic_list, scan_list = [], [], []

    for spectrum in run:
        if spectrum.ms_level == 1:
            rt = float(spectrum.scan_time[0])  # minutes
            mz_arr = np.array(spectrum.mz, dtype=np.float32)
            int_arr = np.array(spectrum.i, dtype=np.float32)
            tic = float(np.sum(int_arr))

            rt_list.append(rt)
            tic_list.append(tic)
            scan_list.append({
                'rt': rt,
                'mz': mz_arr,
                'intensity': int_arr
            })

    mz_all = np.concatenate([s['mz'] for s in scan_list if len(s['mz']) > 0])

    sample_name = Path(mzml_path).stem

    return {
        'rt': np.array(rt_list, dtype=np.float64),
        'tic': np.array(tic_list, dtype=np.float64),
        'scan_list': scan_list,
        'mz_range': (float(mz_all.min()), float(mz_all.max())),
        'sample_name': sample_name,
        'n_scans': len(scan_list)
    }


def load_sample(path: str) -> dict:
    """Load a sample by extension: .RAW read natively (Thermo RawFileReader,
    no conversion), anything else via pymzml. If the native reader's DLLs are
    not available, fall back to RAW->mzML conversion. Same return structure."""
    low = str(path).lower()
    if low.endswith('.qgd'):                     # Shimadzu GCMSsolution
        from qgd_reader import load_qgd_to_matrix
        return load_qgd_to_matrix(path)
    if low.endswith('.raw'):                     # Thermo
        try:
            from raw_reader import load_raw_to_matrix
            return load_raw_to_matrix(path)
        except RuntimeError as e:
            print(f"       native RAW read unavailable ({e}); converting to mzML")
            out_dir = Path(path).parent / "_mzml_tmp"
            mzml = convert_raw_to_mzml(path, str(out_dir))
            return load_mzml_to_matrix(mzml)
    return load_mzml_to_matrix(path)


def extract_eic(scan_list: list, target_mz: float, tolerance: float = 0.5) -> np.ndarray:
    """Extract Extracted Ion Chromatogram for a target m/z.

    Args:
        scan_list: List of scan dicts from load_mzml_to_matrix
        target_mz: Target m/z value
        tolerance: ±Da window (0.5 default for low-res quadrupole)

    Returns:
        1D numpy array of EIC intensities
    """
    eic = []
    for scan in scan_list:
        mask = np.abs(scan['mz'] - target_mz) <= tolerance
        eic.append(float(scan['intensity'][mask].sum()) if mask.any() else 0.0)
    return np.array(eic, dtype=np.float64)
