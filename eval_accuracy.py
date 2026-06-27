#!/usr/bin/env python3
"""Evaluate spectral matching accuracy at user-specified RTs."""
import os, json, openpyxl, re, time, numpy as np

t0 = time.time()
BASE = os.environ['HOME']
txt_path = os.path.join(BASE, 'Desktop', '验证', 'xt6.26.txt')
xlsx_path = os.path.join(BASE, 'Desktop', '验证', 'xt_filled.xlsx')
out_path = os.path.join(BASE, 'Desktop', '验证', 'accuracy_results.json')

import sys
sys.path.insert(0, os.path.join(BASE, 'Desktop', '桌面', 'gcms-pipeline'))
from gcms_pipeline import parse_thermo_txt, load_libraries, pre_search_candidates, \
    nist_similarity, classify_peak, CONTAMINANTS, AQUEOUS_IMPOSSIBLE, Config, \
    load_ri_database, calc_ri_from_rt, compute_ri_score, get_literature_ri, DEFAULT_ALKANE_RTS

# Load library
print('Loading library...')
library, bp_index, ion_masks = load_libraries()

# Load RI database
ri_db = load_ri_database()
ri_column = 'DB-WAX'
ri_tolerance = 50
alkane_rts = DEFAULT_ALKANE_RTS.get(ri_column)
print(f'RI: {ri_column} ±{ri_tolerance}, {len(alkane_rts)} alkane points')
ri_weight = 0.15

print('Parsing TXT...')
rts, tics, spectra = parse_thermo_txt(txt_path)
print(f'{len(rts):,d} scans ({time.time()-t0:.1f}s)')

# Load ground truth
wb = openpyxl.load_workbook(xlsx_path, data_only=True)
ws = wb['Sheet1']
ground_truth = []
for r in range(6, ws.max_row + 1):
    rt = ws.cell(row=r, column=1).value
    name = ws.cell(row=r, column=9).value or ws.cell(row=r, column=10).value or ''
    if rt and name:
        ground_truth.append({'rt': float(rt), 'name': str(name).strip()})

print(f'Ground truth: {len(ground_truth)} compounds')

# Name matching
def is_match(gt_name, pr_name):
    def norm(s):
        return re.sub(r'[^a-z0-9]', '', s.lower().strip())
    gt_n, pr_n = norm(gt_name), norm(pr_name)

    if gt_n == pr_n: return True
    if len(gt_n) > 8 and (gt_n in pr_n or pr_n in gt_n): return True

    # Check key fragments
    gt_frags = set()
    for i in range(len(gt_n)-4):
        frag = gt_n[i:i+5]
        if len(frag) >= 5: gt_frags.add(frag)
    pr_frags = set()
    for i in range(len(pr_n)-4):
        frag = pr_n[i:i+5]
        if len(frag) >= 5: pr_frags.add(frag)
    common = gt_frags & pr_frags
    if len(common) >= 3: return True

    # Known synonym pairs
    known = [
        ('trimethylpyridine', 'pyridine'),
        ('dimethylsilanediol', 'silanediol'),
        ('betaionone', 'ionone'),
        ('trimethylpentanediolmonoisobutyrate', 'trimethylpentyl'),
        ('trimethylpentanedioldiisobutyrate', 'trimethylpentyl'),
        ('nonanal', 'nonanal'),
        ('octanal', 'octanal'),
        ('benzaldehyde', 'benzaldehyde'),
        ('pentylfuran', 'pentylfuran'),
        ('octen3ol', 'octenol'),
        ('octen3one', 'octenone'),
        ('octanedione', 'octanedione'),
        ('heptenal', 'heptenal'),
        ('octenal', 'octenal'),
        ('nonenal', 'nonenal'),
        ('decadienal', 'decadienal'),
        ('farnesyl', 'farnesyl'),
        ('hexadecene', 'hexadecene'),
        ('hexadecenal', 'hexadecenal'),
        ('cyclocitral', 'cyclocitral'),
        ('ditertbutylphenol', 'ditertbutyl'),
        ('butylatedhydroxytoluene', 'butylated'),
        ('octaethylene', 'octaethylene'),
        ('butoxyethoxy', 'butoxy'),
        ('ethylhexanol', 'ethylhexanol'),
        ('dimethylpyrrole', 'pyrrole'),
        ('linalyl', 'linalyl'),
        ('isonone', 'ionone'),
    ]
    for a, b in known:
        if a in gt_n and b in pr_n: return True
        if b in gt_n and a in pr_n: return True
    return False

# Test matching at each user RT
results = []
for gt in ground_truth:
    target_rt = gt['rt']
    # Find closest scan
    idx = np.argmin(np.abs(rts - target_rt))
    rt_diff = abs(rts[idx] - target_rt)

    if rt_diff > 0.2:  # Too far from any scan
        results.append({
            'gt': gt, 'status': 'no_scan',
            'rt_actual': round(rts[idx], 3), 'rt_diff': round(rt_diff, 3)
        })
        continue

    spec = spectra[idx]
    actual_rt = round(rts[idx], 3)
    clean = [(int(round(m)), int(i)) for m, i in spec if int(round(m)) > 25 and i > 100]

    if len(clean) < 5:
        results.append({
            'gt': gt, 'status': 'low_peaks',
            'rt_actual': actual_rt, 'n_peaks': len(clean)
        })
        continue

    tier, peak_label, keep = classify_peak(clean[:15])
    if not keep:
        results.append({
            'gt': gt, 'status': 'L3_rejected',
            'rt_actual': actual_rt, 'tier': tier, 'label': peak_label
        })
        continue

    candidates = pre_search_candidates(clean, bp_index, ion_masks, library)
    if not candidates:
        results.append({
            'gt': gt, 'status': 'no_candidates',
            'rt_actual': actual_rt
        })
        continue

    # Full matching
    threshold = Config.mf_medium if tier == 'L1' else max(Config.threshold, Config.mf_high)
    matches = []

    for idx_c in candidates[:500]:
        comp = library[idx_c]
        mf, rmf, n_sh, n_u, n_l = nist_similarity(clean, comp['peaks'])
        si = max(mf, rmf)
        if si >= threshold:
            matches.append((mf, rmf, si, n_sh, comp['name'], comp.get('cas', '')))

    if not matches:
        results.append({
            'gt': gt, 'status': 'no_match',
            'rt_actual': actual_rt, 'n_candidates': len(candidates)
        })
        continue

    # Filter aqueous-impossible
    matches = [(mf,rmf,si,n,name,cas) for mf,rmf,si,n,name,cas in matches
               if not (name in AQUEOUS_IMPOSSIBLE and max(mf,rmf) < 950)]

    # Phthalate check
    has_149 = any(int(round(m)) == 149 for m, i in clean[:20])
    has_57 = any(int(round(m)) == 57 for m, i in clean[:20])
    if has_149 and has_57:
        results.append({
            'gt': gt, 'status': 'phthalate_blocked',
            'rt_actual': actual_rt
        })
        continue

    matches.sort(key=lambda x: -x[2])

    # ---- RI re-ranking ----
    ri_measured = calc_ri_from_rt(actual_rt, alkane_rts, ri_column)

    reranked = []
    for mf_i, rmf_i, si_i, n_sh_i, name_i, cas_i in matches:
        ri_lit = get_literature_ri(name_i, ri_db, ri_column, cas=cas_i)
        ri_score = compute_ri_score(ri_measured, ri_lit, ri_tolerance)
        # Hard filter: reject if RI deviation > 2x tolerance
        if ri_lit is not None and ri_measured is not None:
            if abs(ri_measured - ri_lit) > 2 * ri_tolerance:
                continue
        combined_si = int(round(si_i * (1 - ri_weight) + ri_score * 999 * ri_weight))
        reranked.append((combined_si, mf_i, rmf_i, si_i, n_sh_i, name_i, cas_i, ri_measured, ri_lit, ri_score))

    if reranked:
        reranked.sort(key=lambda x: -x[0])
    else:
        reranked = [(m[2],) + m for m in matches[:10]]  # fallback without RI
        reranked.sort(key=lambda x: -x[0])

    # Check top 1 and top 5
    top1 = reranked[0]
    is_top1_correct = is_match(gt['name'], top1[5])  # name is now index 5

    top5_names = [m[5] for m in reranked[:5]]
    is_top5_correct = any(is_match(gt['name'], name) for name in top5_names)

    results.append({
        'gt': gt,
        'rt_actual': actual_rt,
        'top1_name': top1[5], 'top1_cas': top1[6],
        'top1_mf': top1[1], 'top1_rmf': top1[2], 'top1_si': top1[3], 'top1_n_shared': top1[4],
        'top1_ri_measured': top1[7], 'top1_ri_lit': top1[8], 'top1_ri_score': top1[9],
        'top5_names': top5_names,
        'correct_top1': is_top1_correct,
        'correct_top5': is_top5_correct,
        'tier': tier,
        'n_candidates': len(candidates),
        'n_matches': len(matches),
    })

# Summary
correct_1 = sum(1 for r in results if r.get('correct_top1'))
correct_5 = sum(1 for r in results if r.get('correct_top5'))
wrong = sum(1 for r in results if r.get('correct_top1') is False)
no_match = sum(1 for r in results if r.get('status') in ('no_match', 'no_candidates'))
rejected = sum(1 for r in results if r.get('status') in ('L3_rejected', 'phthalate_blocked'))
other = sum(1 for r in results if r.get('status') in ('no_scan', 'low_peaks'))

n_total = len(results)

print(f'\n=== MATCHING ACCURACY ===')
print(f'  Top-1 correct:  {correct_1:2d}')
print(f'  Top-5 correct:  {correct_5:2d}')
print(f'  Top-1 wrong:    {wrong:2d}')
print(f'  No match:       {no_match:2d}')
print(f'  Rejected:       {rejected:2d}')
print(f'  Other:          {other:2d}')
print(f'  Total GT:       {n_total:2d}')
print(f'  Acc (top-1): {correct_1}/{correct_1+wrong} = {correct_1/(correct_1+wrong)*100:.0f}%' if correct_1+wrong > 0 else 'N/A')
print(f'  Acc (top-5): {correct_5}/{correct_1+wrong} = {correct_5/(correct_1+wrong)*100:.0f}%' if correct_1+wrong > 0 else 'N/A')

print(f'\n--- CORRECT (Top-1) ---')
for r in results:
    if r.get('correct_top1'):
        gt = r['gt']
        ri_m = r.get('top1_ri_measured', '?')
        ri_l = r.get('top1_ri_lit', '?')
        ri_str = ' RI={}/{}'.format(ri_m, ri_l) if ri_m != '?' else ''
        gt_name = gt['name']
        gt_rt = gt['rt']
        t_si = r['top1_si']
        t_mf = r['top1_mf']
        t_rmf = r['top1_rmf']
        print('  RT={:5.1f} {:45s} SI={:3d} MF={:3d} RMF={:3d}{}'.format(gt_rt, gt_name, t_si, t_mf, t_rmf, ri_str))

print(f'\n--- WRONG (Top-1) ---')
for r in results:
    if 'correct_top1' in r and not r['correct_top1']:
        gt = r['gt']
        print(f'  RT={gt["rt"]:5.1f} GT: {gt["name"]:45s} -> {r["top1_name"]:45s} SI={r["top1_si"]:3d}')

print(f'\n--- NO MATCH ---')
for r in results:
    if r.get('status') in ('no_match', 'no_candidates'):
        gt = r['gt']
        print(f'  RT={gt["rt"]:5.1f} {gt["name"]:45s} ({r.get("status")})')

print(f'\n--- REJECTED ---')
for r in results:
    if r.get('status') in ('L3_rejected', 'phthalate_blocked'):
        gt = r['gt']
        print(f'  RT={gt["rt"]:5.1f} {gt["name"]:45s} ({r.get("label",r.get("status",""))})')

# Save detailed results
serializable = []
for r in results:
    sr = {}
    for k, v in r.items():
        if isinstance(v, dict):
            sr[k] = dict(v)
        elif isinstance(v, (np.integer,)): sr[k] = int(v)
        elif isinstance(v, (np.floating,)): sr[k] = float(v)
        elif isinstance(v, list): sr[k] = [str(x) for x in v]
        else: sr[k] = v
    serializable.append(sr)

with open(out_path, 'w', encoding='utf-8') as f:
    json.dump(serializable, f, ensure_ascii=False, indent=2, default=str)

print(f'\nSaved to {out_path}')
print(f'Total time: {time.time()-t0:.1f}s')
