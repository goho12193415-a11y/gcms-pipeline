#!/usr/bin/env python3
"""Full pipeline vs user manual processing comparison."""
import os, sys, json, re, time, numpy as np, openpyxl

t0 = time.time()
BASE = os.environ['HOME']
TXT = os.path.join(BASE, 'Desktop', '验证', 'xt6.26.txt')
XLSX = os.path.join(BASE, 'Desktop', '验证', 'xt_filled.xlsx')

sys.path.insert(0, os.path.join(BASE, 'Desktop', '桌面', 'gcms-pipeline'))
from gcms_pipeline import (parse_thermo_txt, detect_peaks, load_libraries,
    pre_search_candidates, nist_similarity, classify_peak,
    CONTAMINANTS, AQUEOUS_IMPOSSIBLE, Config, load_ri_database,
    calc_ri_from_rt, get_literature_ri, compute_ri_score, DEFAULT_ALKANE_RTS)

print('=' * 65)
print('  PIPELINE v4.1 vs USER MANUAL PROCESSING — FULL COMPARISON')
print('=' * 65)

# ===== LOAD =====
print('\n[1] Loading...')
library, bp_index, ion_masks = load_libraries()
ri_db = load_ri_database()
alkane_rts = DEFAULT_ALKANE_RTS.get('DB-WAX')
rts, tics, spectra = parse_thermo_txt(TXT)

# ===== USER DATA =====
wb = openpyxl.load_workbook(XLSX, data_only=True)
ws = wb['Sheet1']
user_peaks = []
for r in range(6, ws.max_row + 1):
    rt = ws.cell(row=r, column=1).value
    name = ws.cell(row=r, column=9).value or ws.cell(row=r, column=10).value or ''
    area = ws.cell(row=r, column=4).value or 0
    area_pct = ws.cell(row=r, column=5).value or 0
    cas = ws.cell(row=r, column=8).value or ''
    if rt and str(name).strip():
        user_peaks.append({
            'rt': float(rt), 'name': str(name).strip(), 'cas': str(cas).strip(),
            'area': float(area), 'area_pct': float(area_pct),
        })

# ===== PIPELINE PEAK DETECTION =====
peaks = detect_peaks(rts, tics, spectra)

# ===== PIPELINE IDENTIFICATION =====
pipeline_results = []
for p in peaks:
    raw = p['spectrum']
    clean = [(int(round(m)), int(i)) for m, i in raw if int(round(m)) > 25 and i > 100]
    if len(clean) < 5: continue

    tier, label, keep = classify_peak(clean[:15])
    if not keep: continue

    cand = pre_search_candidates(clean, bp_index, ion_masks, library)
    if not cand: continue

    threshold = Config.mf_medium if tier == 'L1' else max(Config.threshold, Config.mf_high)
    matches = []
    for ci in cand[:500]:
        comp = library[ci]
        mf, rmf, n_sh, n_u, n_l = nist_similarity(clean, comp['peaks'])
        si = max(mf, rmf)
        if si >= threshold:
            matches.append((si, mf, rmf, n_sh, comp))

    if not matches: continue

    matches = [(s,mf,rmf,n,c) for s,mf,rmf,n,c in matches
               if not (c['name'] in AQUEOUS_IMPOSSIBLE and max(mf,rmf)<950)]
    has_149 = any(int(round(m)) == 149 for m, i in clean[:20])
    has_57 = any(int(round(m)) == 57 for m, i in clean[:20])
    if has_149 and has_57: continue

    matches.sort(key=lambda x: -x[0])
    si, mf, rmf, n_sh, comp = matches[0]
    if comp['name'] in CONTAMINANTS: continue

    ri_meas = calc_ri_from_rt(p['rt'], alkane_rts, 'DB-WAX')
    ri_lit = get_literature_ri(comp['name'], ri_db, 'DB-WAX', cas=comp.get('cas',''))

    pipeline_results.append({
        'rt': round(float(p['rt']), 2), 'name': comp['name'], 'cas': comp.get('cas',''),
        'si': si, 'mf': mf, 'rmf': rmf,
        'area': round(float(p['area']), 0), 'tic': round(float(p['tic']), 0),
        'tier': tier, 'ri_meas': ri_meas, 'ri_lit': ri_lit,
    })

# Deduplicate
pipeline_results.sort(key=lambda x: -x['si'])
seen = set()
pipe_final = []
for r in pipeline_results:
    key = (round(r['rt'], 1), r['name'])
    if key not in seen:
        seen.add(key)
        pipe_final.append(r)
pipe_final.sort(key=lambda x: x['rt'])

# ===== COMPARISON =====
def norm(s):
    return re.sub(r'[^a-z0-9]', '', s.lower().strip())

def names_match(a, b):
    an, bn = norm(a), norm(b)
    if an == bn: return True
    if len(an) > 8 and (an in bn or bn in an): return True
    frags_a = {an[i:i+5] for i in range(len(an)-4) if len(an[i:i+5])>=5}
    frags_b = {bn[i:i+5] for i in range(len(bn)-4) if len(bn[i:i+5])>=5}
    if len(frags_a & frags_b) >= 3: return True
    return False

def cas_match(a, b):
    a = re.sub(r'[^0-9-]', '', str(a))
    b = re.sub(r'[^0-9-]', '', str(b))
    return a == b and len(a) > 4

# Match pipe to user by RT (within 0.3 min)
comparisons = []
for up in user_peaks:
    best = None
    best_dist = 999
    for pr in pipe_final:
        d = abs(pr['rt'] - up['rt'])
        if d < best_dist:
            best_dist = d
            best = pr

    status = 'no_peak'
    if best and best_dist <= 0.5:
        if names_match(up['name'], best['name']) or cas_match(up.get('cas',''), best.get('cas','')):
            status = 'agree'
        else:
            status = 'disagree'
    elif best:
        status = 'no_nearby_peak'

    comparisons.append({
        'user_rt': up['rt'], 'user_name': up['name'], 'user_cas': up.get('cas',''),
        'user_area_pct': up.get('area_pct', 0),
        'pipe_rt': best['rt'] if best else None,
        'pipe_name': best['name'] if best else '',
        'pipe_si': best['si'] if best else 0,
        'pipe_mf': best['mf'] if best else 0,
        'pipe_rmf': best['rmf'] if best else 0,
        'pipe_ri_meas': best.get('ri_meas') if best else None,
        'pipe_ri_lit': best.get('ri_lit') if best else None,
        'status': status, 'rt_diff': round(best_dist, 2),
    })

# Find pipe peaks the user missed
user_rts = set(round(u['rt'], 1) for u in user_peaks)
pipe_only = [p for p in pipe_final if all(abs(p['rt']-ur) > 0.5 for ur in user_rts)]

# ===== REPORT =====
agree = [c for c in comparisons if c['status'] == 'agree']
disagree = [c for c in comparisons if c['status'] == 'disagree']
no_peak = [c for c in comparisons if c['status'] == 'no_peak']
no_nearby = [c for c in comparisons if c['status'] == 'no_nearby_peak']

print(f'\n{"="*65}')
print(f'  COMPARISON SUMMARY')
print(f'{"="*65}')
print(f'  User identified peaks:     {len(user_peaks):3d}')
print(f'  Pipeline detected peaks:   {len(peaks):3d}')
print(f'  Pipeline identified:       {len(pipe_final):3d}')
print(f'  Pipeline-only peaks:       {len(pipe_only):3d}  (pipeline found, user missed)')
print(f'')
print(f'  AGREEMENT:                 {len(agree):3d}  (same compound at same RT)')
print(f'  DISAGREEMENT:              {len(disagree):3d}  (different ID at same RT)')
print(f'  PIPELINE MISSED:           {len(no_peak):3d}  (no pipeline peak near user RT)')
print(f'  NO NEARBY:                 {len(no_nearby):3d}')

print(f'\n{"="*65}')
print(f'  AGREEMENT ({len(agree)} peaks)')
print(f'{"="*65}')
for c in sorted(agree, key=lambda x: x['user_rt']):
    ri_str = f'  RI={c["pipe_ri_meas"]}/{c["pipe_ri_lit"]}' if c.get('pipe_ri_lit') else ''
    print(f'  RT={c["user_rt"]:5.1f}  {c["user_name"]:45s}  SI={c["pipe_si"]:3d} MF={c["pipe_mf"]:3d} RMF={c["pipe_rmf"]:3d}{ri_str}')

print(f'\n{"="*65}')
print(f'  DISAGREEMENT ({len(disagree)} peaks)')
print(f'{"="*65}')
for c in sorted(disagree, key=lambda x: x['user_rt']):
    print(f'  RT={c["user_rt"]:5.1f}')
    print(f'    User:   {c["user_name"]}')
    print(f'    Pipe:   {c["pipe_name"]}  (SI={c["pipe_si"]:3d} MF={c["pipe_mf"]:3d} RMF={c["pipe_rmf"]:3d})')

print(f'\n{"="*65}')
print(f'  PIPELINE MISSED ({len(no_peak)} peaks — user found, pipeline no peak)')
print(f'{"="*65}')
for c in sorted(no_peak, key=lambda x: x['user_rt']):
    print(f'  RT={c["user_rt"]:5.1f}  {c["user_name"]:45s}  area={c["user_area_pct"]:.1f}%')

print(f'\n{"="*65}')
print(f'  PIPELINE-ONLY ({len(pipe_only)} peaks — user might have missed)')
print(f'{"="*65}')
for p in sorted(pipe_only, key=lambda x: x['rt'])[:20]:
    ri_str = f'  RI={p.get("ri_meas")}/{p.get("ri_lit")}' if p.get('ri_lit') else ''
    print(f'  RT={p["rt"]:5.1f}  {p["name"]:45s}  SI={p["si"]:3d} MF={p["mf"]:3d} RMF={p["rmf"]:3d}{ri_str} [{p["tier"]}]')
if len(pipe_only) > 20:
    print(f'  ... and {len(pipe_only)-20} more')

elapsed = time.time() - t0
print(f'\n{"="*65}')
print(f'  Elapsed: {elapsed:.1f}s')
print(f'{"="*65}')

# Save comparison
out = os.path.join(BASE, 'Desktop', '验证', 'comparison_report.json')
with open(out, 'w', encoding='utf-8') as f:
    json.dump({
        'agree': len(agree), 'disagree': len(disagree),
        'pipeline_missed': len(no_peak), 'pipeline_only': len(pipe_only),
        'details': comparisons,
        'pipeline_only': [{'rt': p['rt'], 'name': p['name'], 'si': p['si']} for p in pipe_only]
    }, f, ensure_ascii=False, indent=2, default=str)
print(f'Saved to comparison_report.json')
