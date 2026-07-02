#!/usr/bin/env python3
"""
detox_publication_ready.py

Create a publication-ready plot and recommendations from the detox glucose sweep CSV.
"""
from pathlib import Path
import sys

OUT_DIR = Path("micromol_results") / "optimisation" / "detox_glucose_sweep"
CSV_IN = OUT_DIR / "detox_glucose_sweep_summary.csv"
PLOT_OUT = OUT_DIR / "detox_summary_plot.png"
RECS_OUT = OUT_DIR / "detox_recommendations.csv"
TXT_OUT = OUT_DIR / "detox_summary.txt"

if not CSV_IN.exists():
    print(f"Input CSV not found: {CSV_IN}")
    sys.exit(1)

try:
    import pandas as pd
    import matplotlib.pyplot as plt
except Exception as e:
    print("Missing dependency: pandas and matplotlib are required.")
    raise

df = pd.read_csv(CSV_IN)

if 'glucose_mM' not in df.columns:
    raise ValueError("CSV missing 'glucose_mM' column")
if 'peak_mM' not in df.columns:
    # try fallback
    if 'max_mM' in df.columns:
        df = df.rename(columns={'max_mM': 'peak_mM'})
    else:
        raise ValueError("CSV missing 'peak_mM' column")

df = df.sort_values('glucose_mM')

# Plot
plt.figure(figsize=(7,4.5))
plt.plot(df['glucose_mM'], df['peak_mM'], marker='o', linestyle='-', color='#2E8B57', linewidth=2)
plt.scatter(df['glucose_mM'], df['peak_mM'], s=60, color='#2E8B57')
plt.xlabel('Glucose (mM)')
plt.ylabel('Peak detoxified LPS (mM)')
plt.title('Detoxified LPS production vs glucose')
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(PLOT_OUT, dpi=300)
print(f"Saved plot: {PLOT_OUT}")

# Recommendations CSV
rec_df = df[['glucose_mM','peak_mM']].copy()
rec_df = rec_df.rename(columns={'peak_mM':'detoxified_LPS_peak_mM'})
rec_df['recommended_harvest_time'] = 'peak early (~0.2 h simulated); validate 2-4 h in lab'
rec_df.to_csv(RECS_OUT, index=False)
print(f"Saved recommendations: {RECS_OUT}")

# Plain English summary
best_row = rec_df.loc[rec_df['detoxified_LPS_peak_mM'].idxmax()]
summary = f"DETOXIFIED LPS PRODUCTION SUMMARY\n\nBest simulated condition:\n- Glucose: {best_row['glucose_mM']} mM\n- Peak detoxified LPS: {best_row['detoxified_LPS_peak_mM']:.3f} mM\n- Recommended harvest: {best_row['recommended_harvest_time']}\n\nNotes:\n- Detoxified LPS predicted non-toxic and retains immunogenicity in silico.\n- Use MPLA and wild-type LPS as controls in cellular assays.\n\nFiles generated:\n- {PLOT_OUT}\n- {RECS_OUT}\n- {TXT_OUT}\n"

with open(TXT_OUT, 'w') as f:
    f.write(summary)

print(f"Saved summary: {TXT_OUT}")
