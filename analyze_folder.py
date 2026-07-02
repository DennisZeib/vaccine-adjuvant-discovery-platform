"""
Analyze a run folder (micromol_results/*): extract top candidates, create analysis folder, generate time-series and network.
Usage: python analyze_folder.py <run_folder>
"""
import sys
from pathlib import Path
import json
import csv
import matplotlib.pyplot as plt
import networkx as nx
from micromol_core import create_default_microbes, create_metabolite_registry, CombinationSimulator

if len(sys.argv) < 2:
    print('Usage: python analyze_folder.py <run_folder>')
    sys.exit(1)

RUN = Path(sys.argv[1])
if not RUN.exists():
    print('Run folder not found:', RUN)
    sys.exit(1)

DETAILED = RUN / 'micromol_detailed.json'
if not DETAILED.exists():
    print('Detailed JSON not found in run folder:', DETAILED)
    sys.exit(1)

with open(DETAILED, 'r') as f:
    detailed = json.load(f)

ANALYSIS = RUN / 'analysis'
ANALYSIS.mkdir(parents=True, exist_ok=True)

# Build top summary
top_csv = ANALYSIS / 'top_candidates_summary.csv'
rows = []
for entry in detailed:
    combo = entry.get('combination')
    final_concs = entry.get('final_concentrations', {})
    for s in entry.get('scores', []):
        rows.append({'combination': combo, 'molecule_id': s.get('molecule_id'), 'score': s.get('overall_vaccine_score'), 'final_concentration': final_concs.get(s.get('molecule_id'), 0)})

# write summary
with open(top_csv, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=['combination','molecule_id','score','final_concentration'])
    writer.writeheader()
    writer.writerows(rows)

print('Wrote summary to', top_csv)

# Rank and network
import pandas as pd
df = pd.read_csv(top_csv)
df['score'] = pd.to_numeric(df['score'], errors='coerce').fillna(0)
df['final_concentration'] = pd.to_numeric(df['final_concentration'], errors='coerce').fillna(0)
df['rank_metric'] = df['score'] * df['final_concentration']
df_sorted = df.sort_values('rank_metric', ascending=False)
df_sorted.to_csv(ANALYSIS / 'top_candidates_ranked.csv', index=False)
print('Wrote ranked CSV')

# Create network image
G = nx.Graph()
for _, row in df_sorted.iterrows():
    combo = row['combination']
    metric = float(row['rank_metric'])
    parts = [p.strip() for p in combo.split('+')]
    if len(parts) != 2:
        continue
    a,b = parts
    G.add_node(a)
    G.add_node(b)
    if G.has_edge(a,b):
        if metric > G[a][b]['weight']:
            G[a][b]['weight'] = metric
    else:
        G.add_edge(a,b, weight=metric)

png = ANALYSIS / 'top_candidates_network.png'
plt.figure(figsize=(8,8))
pos = nx.spring_layout(G, seed=42)
weights = [G[u][v]['weight'] for u,v in G.edges()]
maxw = max(weights) if weights else 1
widths = [max(0.5, 6*(w/maxw)) for w in weights]
nx.draw_networkx_nodes(G, pos, node_size=800, node_color='#88ccee')
nx.draw_networkx_labels(G, pos)
nx.draw_networkx_edges(G, pos, width=widths)
plt.title('Top candidate combinations')
plt.axis('off')
plt.tight_layout()
plt.savefig(png)
print('Wrote network to', png)

# Re-run time series for top combos (top 10 unique combos)
microbes = create_default_microbes()
metabolites = create_metabolite_registry()
sim = CombinationSimulator(microbes, metabolites)
unique_combos = df_sorted['combination'].unique().tolist()[:10]
for combo in unique_combos:
    combo_key = combo
    parts = [p.strip() for p in combo.split('+')]
    parts = [p for p in parts if p]
    print('Re-running', parts)
    res = sim.simulate_coculture(parts, {'glucose':20.0}, duration=24.0)
    ch = res.get('concentration_history')
    if ch:
        times = res.get('time_steps', [])
        keys = sorted(ch.keys())
        out_csv = ANALYSIS / f"{combo.replace(' ','').replace('+','_')}_timeseries.csv"
        with open(out_csv, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['time'] + keys)
            for i,t in enumerate(times):
                writer.writerow([t] + [ch[k][i] for k in keys])
        print(' Wrote timeseries', out_csv)

print('Analysis complete for', RUN)
