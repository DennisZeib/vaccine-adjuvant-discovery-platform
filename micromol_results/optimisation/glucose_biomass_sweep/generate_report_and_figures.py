"""
Generate publication-ready figures and a CSV report (harvest times + top candidates)
from sweep outputs produced by optimise_glucose_and_biomass.py
"""
import json
import csv
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

plt.style.use('ggplot')

BASE = Path('.')
OUTDIR = BASE / 'micromol_results' / 'optimisation' / 'glucose_biomass_sweep'
if not OUTDIR.exists():
    raise SystemExit(f"Output directory not found: {OUTDIR}")

summary_csv = OUTDIR / 'optimisation_summary.csv'
if not summary_csv.exists():
    raise SystemExit(f"Summary CSV not found: {summary_csv}")

df = pd.read_csv(summary_csv)
# Normalize microbe pair display
pairs = df[['microA','microB']].drop_duplicates().values.tolist()

# Plot: LPS vs Glucose for each pair, colored by ratio
for microA, microB in pairs:
    sel = df[(df.microA==microA)&(df.microB==microB)]
    if sel.empty:
        continue
    plt.figure(figsize=(8,6))
    for r, grp in sel.groupby('ratio_microA'):
        plt.plot(grp['glucose_mM'], grp['lps_mM'], marker='o', label=str(r))
    plt.title(f"LPS production: {microA} + {microB}")
    plt.xlabel('Glucose (mM)')
    plt.ylabel('LPS (mM)')
    plt.legend(title=f'{microA} : {microB} ratio')
    plt.tight_layout()
    outpng = OUTDIR / f"{microA}_{microB}_LPS_vs_glucose.png"
    plt.savefig(outpng, dpi=300)
    plt.close()

# Plot: LPS vs Ratio for each glucose level
for microA, microB in pairs:
    sel = df[(df.microA==microA)&(df.microB==microB)]
    if sel.empty:
        continue
    glevels = sorted(sel.glucose_mM.unique())
    for glc in glevels:
        s2 = sel[sel.glucose_mM==glc]
        plt.figure(figsize=(8,6))
        plt.plot(s2['ratio_microA'], s2['lps_mM'], marker='o')
        plt.xscale('log')
        plt.title(f'LPS vs {microA} ratio @ {glc} mM glucose ({microA}+{microB})')
        plt.xlabel(f'{microA} : {microB} biomass ratio')
        plt.ylabel('LPS (mM)')
        plt.tight_layout()
        outpng = OUTDIR / f"{microA}_{microB}_LPS_vs_ratio_glc{int(glc)}.png"
        plt.savefig(outpng, dpi=300)
        plt.close()

# Now analyze time-series from detailed JSONs to recommend harvest times
records = []
for jf in OUTDIR.glob('*.json'):
    try:
        data = json.load(jf.open())
    except Exception:
        continue
    pair = data.get('pair', [])
    if not pair or len(pair) < 2:
        continue
    microA, microB = pair[0], pair[1]
    glucose = data.get('glucose', None)
    ratio = data.get('ratio', {})
    res = data.get('result', {})
    time_steps = res.get('time_steps', [])
    ch = res.get('concentration_history', {})
    lps_ts = ch.get('polysaccharide_lps', [])
    if not lps_ts:
        max_lps = 0.0
        t_max = None
        t_95 = None
    else:
        arr = np.array(lps_ts)
        max_lps = float(arr.max())
        t_idx = int(arr.argmax())
        t_max = float(time_steps[t_idx]) if t_idx < len(time_steps) else None
        # time to reach 95% of max
        thresh = 0.95 * max_lps
        try:
            idx95 = int(np.where(arr >= thresh)[0][0])
            t_95 = float(time_steps[idx95])
        except Exception:
            t_95 = t_max
    records.append({
        'microA': microA,
        'microB': microB,
        'glucose_mM': glucose,
        'ratio_microA': ratio.get(microA, ''),
        'max_lps_mM': max_lps,
        't_max_h': t_max,
        't_95_h': t_95,
        'json_file': jf.name
    })

recdf = pd.DataFrame(records)

if not recdf.empty:
    # Save harvest times
    recdf.to_csv(OUTDIR / 'harvest_times.csv', index=False)
    # Top candidates overall
    top = recdf.sort_values('max_lps_mM', ascending=False).head(20)
    top.to_csv(OUTDIR / 'top_candidates.csv', index=False)

    # Also create a summary plot: heatmap of max LPS for each pair x glucose
    for microA, microB in pairs:
        sub = recdf[(recdf.microA==microA)&(recdf.microB==microB)]
        if sub.empty:
            continue
        pivot = sub.pivot_table(index='ratio_microA', columns='glucose_mM', values='max_lps_mM')
        plt.figure(figsize=(8,6))
        im = plt.imshow(pivot.values, aspect='auto', cmap='viridis')
        plt.colorbar(im, label='Max LPS (mM)')
        plt.yticks(range(len(pivot.index)), [str(x) for x in pivot.index])
        plt.xticks(range(len(pivot.columns)), [str(int(x)) for x in pivot.columns])
        for (i, j), val in np.ndenumerate(pivot.values):
            plt.text(j, i, f"{val:.3f}", ha='center', va='center', color='white', fontsize=8)
        plt.title(f'Max LPS heatmap: {microA} + {microB}')
        plt.xlabel('Glucose (mM)')
        plt.ylabel(f'{microA} ratio')
        plt.tight_layout()
        plt.savefig(OUTDIR / f"{microA}_{microB}_maxLPS_heatmap.png", dpi=300)
        plt.close()

print(f"Report and figures generated in {OUTDIR}")
