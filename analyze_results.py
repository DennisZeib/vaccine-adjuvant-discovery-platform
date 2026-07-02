"""
Analyze micromol results: extract top vaccine candidates, re-run their simulations to capture time-series, and save CSVs and plots.
"""
import json
from pathlib import Path
import os
import matplotlib.pyplot as plt
import csv

from micromol_core import create_default_microbes, create_metabolite_registry, CombinationSimulator

RESULT_DIR = Path("micromol_results/long_run")
ANALYSIS_DIR = RESULT_DIR / "analysis"
ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)

# Load detailed results
with open(RESULT_DIR / "micromol_detailed.json", "r") as f:
    detailed = json.load(f)

# Find top molecules (overall_vaccine_score >= 70)
top_threshold = 70
candidates = []  # tuples (combo, molecule_id, score, final_conc)
for entry in detailed:
    combo = entry.get("combination")
    final_concs = entry.get("final_concentrations", {})
    for score in entry.get("scores", []):
        if score.get("overall_vaccine_score", 0) >= top_threshold:
            mid = score.get("molecule_id")
            candidates.append((combo, mid, score.get("overall_vaccine_score"), final_concs.get(mid, 0)))

# Write summary CSV
summary_csv = ANALYSIS_DIR / "top_candidates_summary.csv"
with open(summary_csv, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["combination", "molecule_id", "score", "final_concentration"])
    writer.writeheader()
    for combo, mid, sc, conc in candidates:
        writer.writerow({"combination": combo, "molecule_id": mid, "score": sc, "final_concentration": conc})

print(f"Wrote summary to {summary_csv}")

# Re-run simulations for each unique combination in candidates to capture time-series
microbes_db = create_default_microbes()
metabolites = create_metabolite_registry()
sim = CombinationSimulator(microbes_db, metabolites)

unique_combos = sorted({c for c,_,_,_ in candidates})

for combo in unique_combos:
    microbe_ids = [m.strip() for m in combo.split("+")]
    microbe_ids = [m for m in microbe_ids if m]
    microbe_ids = [m for m in microbe_ids]
    print(f"Re-running: {combo} -> {microbe_ids}")

    # run with small dt to capture time series
    network_result = sim.simulate_coculture(microbe_ids, {"glucose": 20.0}, duration=24.0)

    # If concentration_history exists in network_result, save it; otherwise, save final concentrations
    ch = network_result.get("concentration_history")
    if ch:
        # write CSV with time on first column
        times = network_result.get("time_steps", [])
        keys = sorted(ch.keys())
        out_csv = ANALYSIS_DIR / f"{combo.replace(' ', '').replace('+','_')}_timeseries.csv"
        with open(out_csv, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["time"] + keys)
            for i, t in enumerate(times):
                row = [t] + [ch[k][i] if i < len(ch[k]) else "" for k in keys]
                writer.writerow(row)
        print(f"  Wrote time series to {out_csv}")

        # Plot top candidate molecules in this combo
        # find molecules of interest for this combo (those in scores)
        mols = [s.get("molecule_id") for s in (next((e for e in detailed if e.get("combination")==combo),{})).get("scores",[])]
        # select top mols with score >= top_threshold
        top_mols = [s.get("molecule_id") for s in (next((e for e in detailed if e.get("combination")==combo),{})).get("scores",[]) if s.get("overall_vaccine_score",0) >= top_threshold]
        plot_mols = top_mols if top_mols else mols[:3]

        plt.figure(figsize=(8,4))
        for mid in plot_mols:
            if mid in ch:
                plt.plot(times, ch[mid], label=mid)
        plt.xlabel('Time (h)')
        plt.ylabel('Concentration (mM)')
        plt.title(f'Time series - {combo}')
        plt.legend()
        plt.tight_layout()
        out_png = ANALYSIS_DIR / f"{combo.replace(' ', '').replace('+','_')}_timeseries.png"
        plt.savefig(out_png)
        plt.close()
        print(f"  Wrote plot to {out_png}")
    else:
        # Save final concentrations instead
        out_csv = ANALYSIS_DIR / f"{combo.replace(' ', '').replace('+','_')}_final_concentrations.csv"
        with open(out_csv, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["metabolite", "concentration"])
            for k,v in network_result.get("final_concentrations", {}).items():
                writer.writerow([k, v])
        print(f"  Wrote final concentrations to {out_csv}")

print("Analysis complete.")
