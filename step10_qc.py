"""
Step 10: Batch QC Audit
========================
Checks RT drift, area reproducibility, and response factor stability
across multiple samples.
"""
import numpy as np
import pandas as pd


def batch_qc_check(
    all_results: dict,
    qc_sample_pattern: str = "QC",
    rt_drift_threshold: float = 0.05,
    area_cv_threshold: float = 0.15
) -> dict:
    """Run batch-wide QC audit.

    Checks:
      1. RT drift: same compound RT variability across samples
      2. Area reproducibility: CV% across QC replicates
      3. Response factor stability (placeholder)

    Args:
        all_results: {sample_name: DataFrame} from compile_results
        qc_sample_pattern: Substring to identify QC samples
        rt_drift_threshold: Max allowed RT std dev (min)
        area_cv_threshold: Max allowed area CV (%)

    Returns:
        qc_report dict with status, warnings, and errors
    """
    qc_report = {
        'status': 'PASS',
        'warnings': [],
        'errors': [],
        'rt_drift': {},
        'area_cv': {}
    }

    if len(all_results) < 2:
        qc_report['status'] = 'SINGLE_SAMPLE'
        return qc_report

    # Concatenate all results
    all_df = pd.concat(all_results.values(), ignore_index=True)

    # 1. RT drift check (exclude Unknown)
    known = all_df[all_df['Compound_Name'] != 'Unknown']
    for compound, group in known.groupby('Compound_Name'):
        if len(group) < 2:
            continue
        rt_std = group['RT_min'].std()
        if rt_std > rt_drift_threshold:
            qc_report['warnings'].append(
                f"RT drift: {compound} σ={rt_std:.3f} min "
                f"(threshold: {rt_drift_threshold})"
            )
            qc_report['rt_drift'][compound] = round(float(rt_std), 4)

    # 2. QC sample area reproducibility
    qc_samples = {k: v for k, v in all_results.items()
                  if qc_sample_pattern in k}
    if qc_samples:
        qc_all = pd.concat(qc_samples.values(), ignore_index=True)
        for compound, group in qc_all.groupby('Compound_Name'):
            if len(group) < 2 or group['Area'].mean() == 0:
                continue
            cv = group['Area'].std() / group['Area'].mean()
            if cv > area_cv_threshold:
                qc_report['warnings'].append(
                    f"Area CV: {compound} CV={cv*100:.1f}% "
                    f"(threshold: {area_cv_threshold*100:.0f}%)"
                )
                qc_report['area_cv'][compound] = round(float(cv), 4)

    # Determine final status
    if qc_report['errors']:
        qc_report['status'] = 'FAIL'
    elif qc_report['warnings']:
        qc_report['status'] = 'WARN'

    return qc_report
