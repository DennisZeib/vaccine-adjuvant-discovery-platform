#!/usr/bin/env python3
"""Adjuvant sweep: find culture conditions that maximise safe adjuvant yield while minimising toxic byproducts.

Uses `micromol_core.CombinationSimulator` and default microbes/metabolites.
Outputs:
 - bgc_analysis/reports/adjuvant_sweep_results.csv
 - bgc_analysis/reports/top5_adjuvant_sweep.png

Run: python bgc_analysis/adjuvant_sweep.py
"""
import os
from pathlib import Path
import csv
import math

try:
    import pandas as pd
except Exception:
    pd = None

try:
    import matplotlib.pyplot as plt
except Exception:
    plt = None

from micromol_core import create_default_microbes, create_metabolite_registry, CombinationSimulator

# --- Parameters -------------------------------------------------------------
OUTPUT_DIR = Path(__file__).resolve().parent / 'reports'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
CSV_OUT = OUTPUT_DIR / 'adjuvant_sweep_results.csv'
PLOT_OUT = OUTPUT_DIR / 'top5_adjuvant_sweep.png'

TARGET_ADJUVANTS = ['polysaccharide_lps', 'mannan_polysaccharide']
TOXIC_METABOLITES = {
    'ethanol': 1.0,
    'lactate': 0.8,
    'butyrate': 0.6,
}

GLUCOSE_CONCS = [10.0, 20.0, 50.0, 100.0]  # mM
RATIO_VALUES = [1.0, 5.0, 10.0, 20.0]     # salmonella : clostridium (relative biomass)
MICROBES_PAIR = ['salmonella', 'clostridium']
DURATION_HOURS = 48.0

# weighting for overall metric: prefer yield, penalise toxicity
TOXICITY_PENALTY = 0.1

# Create simulator
microbes = create_default_microbes()
metabolites = create_metabolite_registry()
sim = CombinationSimulator(microbes, metabolites)

results = []

print(f"Starting adjuvant sweep for pair: {MICROBES_PAIR} for {DURATION_HOURS} h")
for glucose in GLUCOSE_CONCS:
    for ratio in RATIO_VALUES:
        # biomass ratios: salmonella gets `ratio`, clostridium gets 1.0
        biomass = {MICROBES_PAIR[0]: float(ratio), MICROBES_PAIR[1]: 1.0}
        initial_nutrients = {'glucose': float(glucose)}
        try:
            res = sim.simulate_coculture(MICROBES_PAIR, initial_nutrients,
                                         duration=DURATION_HOURS,
                                         biomass_ratios=biomass)
        except Exception as e:
            print(f"Simulation failed for glucose={glucose}, ratio={ratio}: {e}")
            continue

        final = res.get('final_concentrations', {})
        # compute adjuvant yield (sum over target adjuvants)
        target_yield = 0.0
        for t in TARGET_ADJUVANTS:
            target_yield += float(final.get(t, 0.0))

        # compute toxicity score as weighted sum
        toxicity_score = 0.0
        for tox, weight in TOXIC_METABOLITES.items():
            toxicity_score += float(final.get(tox, 0.0)) * weight

        # overall metric: target_yield penalised by toxicity
        overall_score = target_yield - TOXICITY_PENALTY * toxicity_score

        results.append({
            'glucose_mM': glucose,
            'ratio_salmonella_to_clostridium': f"{int(ratio)}:1",
            'ratio_value': ratio,
            'target_yield_mM': target_yield,
            'toxicity_score': toxicity_score,
            'overall_score': overall_score,
            'final_concentrations': final,
        })
        print(f"glucose={glucose} mM, ratio={int(ratio)}:1 -> yield={target_yield:.4f}, tox={toxicity_score:.4f}, score={overall_score:.4f}")

# Save results to CSV
# flatten final_concentrations for key metabolites we care about
flat_rows = []
for r in results:
    row = {
        'glucose_mM': r['glucose_mM'],
        'ratio_salmonella_to_clostridium': r['ratio_salmonella_to_clostridium'],
        'target_yield_mM': r['target_yield_mM'],
        'toxicity_score': r['toxicity_score'],
        'overall_score': r['overall_score'],
    }
    # include individual target adjuvants
    for t in TARGET_ADJUVANTS:
        row[f'adjuvant_{t}'] = r['final_concentrations'].get(t, 0.0)
    # include toxic metabolites
    for tox in TOXIC_METABOLITES.keys():
        row[f'toxic_{tox}'] = r['final_concentrations'].get(tox, 0.0)
    flat_rows.append(row)

# sort by overall_score desc
flat_rows.sort(key=lambda x: x['overall_score'], reverse=True)

# assign rank
for i, rr in enumerate(flat_rows, start=1):
    rr['overall_rank'] = i

# write CSV
fieldnames = list(flat_rows[0].keys()) if flat_rows else ['glucose_mM','ratio_salmonella_to_clostridium','target_yield_mM','toxicity_score','overall_score','overall_rank']
with open(CSV_OUT, 'w', newline='') as fh:
    writer = csv.DictWriter(fh, fieldnames=fieldnames)
    writer.writeheader()
    for rr in flat_rows:
        writer.writerow(rr)

print(f"Wrote sweep results to {CSV_OUT}")

# Make bar plot of top 5
if plt is not None and len(flat_rows) > 0:
    top5 = flat_rows[:5]
    labels = [f"G{int(r['glucose_mM'])}-R{r['ratio_salmonella_to_clostridium']}" for r in top5]
    scores = [r['overall_score'] for r in top5]
    yields = [r['target_yield_mM'] for r in top5]
    tox = [r['toxicity_score'] for r in top5]

    fig, ax = plt.subplots(figsize=(8,5))
    x = range(len(labels))
    ax.bar(x, yields, color='tab:green', alpha=0.8, label='Adjuvant yield (mM)')
    ax.bar(x, tox, bottom=yields, color='tab:red', alpha=0.6, label='Toxicity (weighted)')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel('mM / weighted units')
    ax.set_title('Top 5 adjuvant-producing conditions (yield + toxicity)')
    ax.legend()
    plt.tight_layout()
    fig.savefig(PLOT_OUT)
    print(f"Saved plot to {PLOT_OUT}")
else:
    print("matplotlib not available; skipped plot")

# If pandas present also print a small table
if pd is not None and len(flat_rows) > 0:
    df = pd.DataFrame(flat_rows)
    print('\nTop 5 conditions:')
    print(df[['glucose_mM','ratio_salmonella_to_clostridium','target_yield_mM','toxicity_score','overall_score']].head(5).to_string(index=False))

print('Done.')
