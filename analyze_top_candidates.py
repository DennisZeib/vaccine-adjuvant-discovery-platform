"""
Analyze top candidates: inspect a timeseries, rank by score * final_concentration, and draw network.
"""
import csv
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx

BASE = Path('micromol_results/long_run/analysis')
BASE.mkdir(parents=True, exist_ok=True)
TOP_CSV = BASE / 'top_candidates_summary.csv'
RANKED_CSV = BASE / 'top_candidates_ranked.csv'
NETWORK_PNG = BASE / 'top_candidates_network.png'

# 1) Inspect a timeseries CSV (ecoli_yeast as example)
ts_file = Path('micromol_results/long_run/analysis/ecoli_yeast_timeseries.csv')
if ts_file.exists():
    df_ts = pd.read_csv(ts_file)
    # show last few rows for LPS and mannan columns if present
    cols_of_interest = [c for c in df_ts.columns if 'polysaccharide_lps' in c or 'mannan_polysaccharide' in c]
    if cols_of_interest:
        sample_summary = df_ts[cols_of_interest].tail(5)
    else:
        sample_summary = None
else:
    df_ts = None

# 2) Rank combinations by score * final_concentration
if TOP_CSV.exists():
    df = pd.read_csv(TOP_CSV)
    # Ensure numeric
    df['score'] = pd.to_numeric(df['score'], errors='coerce').fillna(0)
    df['final_concentration'] = pd.to_numeric(df['final_concentration'], errors='coerce').fillna(0)
    df['rank_metric'] = df['score'] * df['final_concentration']
    df_sorted = df.sort_values('rank_metric', ascending=False)
    df_sorted.to_csv(RANKED_CSV, index=False)
    top10 = df_sorted.head(10)
else:
    top10 = None

# 3) Build network of microbes (edges weighted by rank_metric) using top entries
if TOP_CSV.exists():
    G = nx.Graph()
    for _, row in df_sorted.iterrows():
        combo = row['combination']
        metric = float(row['rank_metric'])
        # parse microbes
        parts = [p.strip() for p in combo.split('+')]
        if len(parts) != 2:
            continue
        a, b = parts
        G.add_node(a)
        G.add_node(b)
        # if multiple edges, keep max
        if G.has_edge(a,b):
            if metric > G[a][b]['weight']:
                G[a][b]['weight'] = metric
        else:
            G.add_edge(a,b, weight=metric)

    # Draw network with edge width scaled
    plt.figure(figsize=(8,8))
    pos = nx.spring_layout(G, seed=42)
    weights = [G[u][v]['weight'] for u,v in G.edges()]
    # scale widths
    maxw = max(weights) if weights else 1
    widths = [max(0.5, 6*(w/maxw)) for w in weights]
    nx.draw_networkx_nodes(G, pos, node_size=800, node_color='#88ccee')
    nx.draw_networkx_labels(G, pos)
    nx.draw_networkx_edges(G, pos, width=widths)
    plt.title('Top candidate combinations (edge width = score * concentration)')
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(NETWORK_PNG)

# Print concise summary
print('Ranked results written to:', RANKED_CSV)
if top10 is not None:
    print('\nTop entries (by score * concentration):')
    print(top10[['combination','molecule_id','score','final_concentration','rank_metric']].to_string(index=False))
if df_ts is not None:
    print('\nTimeseries sample (last rows) for ecoli_yeast:')
    if cols_of_interest:
        print(df_ts[cols_of_interest].tail())
    else:
        print('  No LPS/mannan columns found in timeseries file.')
print('\nNetwork image saved to:', NETWORK_PNG)
