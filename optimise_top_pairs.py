"""
Optimise top candidate pairs across multiple carbon sources.
- Reads top candidates from micromol_results/expanded_run/analysis/top_candidates_ranked.csv
- Selects top N unique combinations
- Runs `run_simulator.py` for each pair and substrate
- Aggregates final concentrations for the target molecule and writes summary CSV + plots
"""
import csv
import subprocess
import sys
from pathlib import Path
import time
import json
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path('micromol_results')
EXPANDED_ANALYSIS = ROOT / 'expanded_run' / 'analysis'
RANKED = EXPANDED_ANALYSIS / 'top_candidates_ranked.csv'
OUTDIR = ROOT / 'optimisation'
OUTDIR.mkdir(parents=True, exist_ok=True)

# Parameters
TOP_N = 5
SUBSTRATES = ['glucose', 'acetate', 'lactate', 'glycerol']
DURATION = 48
CONC = 20.0
SERIAL = True

if not RANKED.exists():
    print('Ranked file not found:', RANKED)
    sys.exit(1)

df = pd.read_csv(RANKED)
# Keep combinations with non-zero rank_metric and sort
if 'rank_metric' in df.columns:
    df_sorted = df.sort_values('rank_metric', ascending=False)
else:
    df_sorted = df.sort_values('score', ascending=False)

# Unique combos ordered
unique_combos = []
for combo in df_sorted['combination'].tolist():
    if combo not in unique_combos:
        unique_combos.append(combo)
    if len(unique_combos) >= TOP_N:
        break

print('Top combos to optimise:', unique_combos)

summary_rows = []

for combo in unique_combos:
    microbe_list = [m.strip() for m in combo.split('+')]
    # Find target metabolite for this combo from ranked df
    row = df[df['combination'] == combo].iloc[0]
    target_met = row['molecule_id']

    for substrate in SUBSTRATES:
        out_dir = OUTDIR / f"{combo.replace(' ', '').replace('+','_')}" / substrate
        out_dir.mkdir(parents=True, exist_ok=True)
        cmd = [sys.executable, 'run_simulator.py', '--microbes'] + microbe_list + [
            '--substrate', substrate,
            '--concentration', str(CONC),
            '--duration', str(DURATION),
            '--serial',
            '--output-dir', str(out_dir)
        ]
        print('Running:', ' '.join(cmd))
        start = time.time()
        subprocess.run(cmd, check=True)
        elapsed = time.time() - start
        print(f'  Completed in {elapsed:.1f}s -> {out_dir}')

        # Read detailed results
        detailed = out_dir / 'micromol_detailed.json'
        final_conc = None
        if detailed.exists():
            with open(detailed, 'r') as f:
                det = json.load(f)
            # det is list; find the entry with this combination
            combo_name = ' + '.join(microbe_list)
            match = None
            for e in det:
                if e.get('combination') == combo_name:
                    match = e
                    break
            if match:
                final_conc = match.get('final_concentrations', {}).get(target_met, 0.0)
            else:
                final_conc = 0.0
        else:
            final_conc = None

        summary_rows.append({
            'combination': combo,
            'microbes': combo_name,
            'target_metabolite': target_met,
            'substrate': substrate,
            'final_concentration': final_conc
        })

# Write summary CSV
summary_csv = OUTDIR / 'optimisation_summary.csv'
with open(summary_csv, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=summary_rows[0].keys())
    writer.writeheader()
    writer.writerows(summary_rows)

print('Wrote summary:', summary_csv)

# Plot results for each combo
for combo in unique_combos:
    dfc = pd.DataFrame([r for r in summary_rows if r['combination']==combo])
    plt.figure()
    plt.bar(dfc['substrate'], [0 if v is None else v for v in dfc['final_concentration']], color='C0')
    plt.xlabel('Substrate')
    plt.ylabel('Final concentration (mM)')
    plt.title(combo + ' - ' + dfc.iloc[0]['target_metabolite'])
    plt.tight_layout()
    outpng = OUTDIR / f"{combo.replace(' ','').replace('+','_')}_substrate_sweep.png"
    plt.savefig(outpng)
    plt.close()
    print('Wrote plot', outpng)

print('Optimisation complete.')
