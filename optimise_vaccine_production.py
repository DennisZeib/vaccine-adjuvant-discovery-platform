#!/usr/bin/env python3
"""
optimise_vaccine_production.py

Sweep glucose concentration to optimise production of polysaccharide_lps
by salmonella + clostridium, while minimising toxic byproducts.

Outputs:
  optimisation_results.csv
  optimisation_plot.png

Usage:
  python optimise_vaccine_production.py
"""

import sys
from pathlib import Path
import csv
import numpy as np

# Ensure project root is in path
_repo_root = Path(__file__).resolve().parent
if str(_repo_root) not in sys.path:
    sys.path.insert(0, str(_repo_root))

from micromol_core import create_default_microbes, create_metabolite_registry, CombinationSimulator

# Optional plotting
try:
    import matplotlib.pyplot as plt
    HAS_PLT = True
except ImportError:
    HAS_PLT = False
    print("matplotlib not installed – plot will be skipped.")

# ----------------------------------------------------------------------
# Parameters
# ----------------------------------------------------------------------
OUTPUT_DIR = Path("bgc_analysis") / "reports"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

CSV_OUT = OUTPUT_DIR / "optimisation_results.csv"
PLOT_OUT = OUTPUT_DIR / "optimisation_plot.png"

MICROBE_A = "salmonella"
MICROBE_B = "clostridium"
TARGET = "polysaccharide_lps"
TOXIC = ["ethanol", "lactate", "acetate"]

GLUCOSE_CONCS = [10, 20, 50, 100]   # mM
SIM_DURATION = 48.0                 # hours

# Biomass ratio pairs (microbe A : microbe B)
BIOMASS_RATIOS = [ (1, 1), (1, 5), (1, 20) ]

# Use top combos from previous substrate sweep. If the optimisation summary
# exists, we'll read combos from it; otherwise fall back to a default list.
OPTIM_SUMMARY = Path("micromol_results") / "optimisation" / "optimisation_summary.csv"

# ----------------------------------------------------------------------
# Setup
# ----------------------------------------------------------------------
microbes = create_default_microbes()
metabolites = create_metabolite_registry()
sim = CombinationSimulator(microbes, metabolites)

# Determine combos to run
if OPTIM_SUMMARY.exists():
    combos = []
    with open(OPTIM_SUMMARY, "r", newline="") as f:
        reader = csv.DictReader(f)
        for r in reader:
            combo = r.get("combination")
            substrate = r.get("substrate")
            # pick combos previously tested (we'll focus on those)
            if combo and combo not in combos and substrate == "glucose":
                combos.append(combo)
    if not combos:
        combos = [f"{MICROBE_A} + {MICROBE_B}"]
else:
    combos = [f"{MICROBE_A} + {MICROBE_B}"]

print("Combos to run:", combos)

# ----------------------------------------------------------------------
# Run glucose sweeps
# ----------------------------------------------------------------------
records = []

# For each combo, sweep glucose concentrations and biomass ratios
for combo in combos:
    combo_safe = combo.replace(" ", "_").replace("+", "_")
    microbes_pair = [m.strip() for m in combo.split("+")]
    if len(microbes_pair) != 2:
        microbes_pair = [MICROBE_A, MICROBE_B]

    combo_out = Path("micromol_results") / "optimisation" / combo_safe
    combo_out.mkdir(parents=True, exist_ok=True)

    for glucose in GLUCOSE_CONCS:
        for (rA, rB) in BIOMASS_RATIOS:
            ratio_str = f"{int(rA)}_{int(rB)}"
            print(f"Running {combo} | glucose={glucose} mM | ratio={ratio_str} ...")
            initial = {"glucose": float(glucose)}
            biomass_map = {microbes_pair[0]: float(rA), microbes_pair[1]: float(rB)}
            try:
                res = sim.simulate_coculture(microbes_pair, initial, duration=SIM_DURATION, biomass_ratios=biomass_map)
            except Exception as e:
                print(f"  Simulation failed: {e}")
                continue

            # Write detailed JSON per-run
            run_dir = combo_out / f"glucose_{int(glucose)}_ratio_{ratio_str}"
            run_dir.mkdir(parents=True, exist_ok=True)
            import json
            with open(run_dir / "micromol_detailed.json", "w") as jf:
                json.dump(res, jf, indent=2)

            # extract timecourse for target and toxicity
            time_steps = res.get("time_steps", [])
            history = res.get("concentration_history", {})
            lps_hist = history.get(TARGET, [])
            peak_val = max(lps_hist) if lps_hist else 0.0
            peak_idx = lps_hist.index(peak_val) if lps_hist else 0
            peak_time = time_steps[peak_idx] if time_steps else 0.0
            final_lps = res.get("final_concentrations", {}).get(TARGET, 0.0)
            tox_max = max((max(history.get(t, [0])) for t in TOXIC), default=0.0)

            records.append({
                "combination": combo,
                "microbes": " + ".join(microbes_pair),
                "glucose_mM": glucose,
                "ratio": f"{int(rA)}:{int(rB)}",
                "peak_lps_mM": peak_val,
                "peak_time_h": peak_time,
                "final_lps_mM": final_lps,
                "toxicity_max_mM": tox_max
            })

            # save a small PNG timecourse per run if plotting available
            if HAS_PLT and lps_hist:
                plt.figure(figsize=(6, 3))
                plt.plot(time_steps, lps_hist, label=f"LPS ({glucose} mM, {rA}:{rB})")
                plt.xlabel("Time (h)")
                plt.ylabel("LPS (mM)")
                plt.title(f"{combo} — glucose {glucose} mM — ratio {rA}:{rB}")
                plt.tight_layout()
                plt.savefig(run_dir / "lps_timecourse.png", dpi=150)
                plt.close()

# ----------------------------------------------------------------------
# Compute toxicity penalty and overall rank (based on peak LPS)
# ----------------------------------------------------------------------
# Normalize toxicity by max observed toxicity_max_mM
max_tox = max((r.get("toxicity_max_mM", 0.0) for r in records), default=0.0)
for r in records:
    if max_tox > 0:
        r["toxicity_penalty"] = r.get("toxicity_max_mM", 0.0) / max_tox
    else:
        r["toxicity_penalty"] = 0.0

# Normalize peak LPS by max and compute overall rank
max_peak = max((r.get("peak_lps_mM", 0.0) for r in records), default=0.0)
for r in records:
    if max_peak > 0:
        norm_peak = r.get("peak_lps_mM", 0.0) / max_peak
        r["overall_rank"] = norm_peak - 0.5 * r.get("toxicity_penalty", 0.0)
    else:
        r["overall_rank"] = float("nan")

# Sort by overall_rank (descending)
records.sort(key=lambda x: x.get("overall_rank", float("-inf")), reverse=True)

# ----------------------------------------------------------------------
# Write CSV
# ----------------------------------------------------------------------
fieldnames = [
    "combination",
    "microbes",
    "glucose_mM",
    "ratio",
    "peak_lps_mM",
    "peak_time_h",
    "final_lps_mM",
    "toxicity_max_mM"
]
# include scoring fields
fieldnames += ["toxicity_penalty", "overall_rank"]

final_out = Path("micromol_results") / "optimisation" / "glucose_biomass_sweep" / "optimisation_summary.csv"
final_out.parent.mkdir(parents=True, exist_ok=True)
with open(final_out, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    for r in records:
        writer.writerow(r)

print(f"Combined results saved to {final_out}")

# ----------------------------------------------------------------------
# Plot top 5 conditions
# ----------------------------------------------------------------------
if HAS_PLT:
    top5 = records[:5]
    labels = [f"{int(r['glucose_mM'])} mM" for r in top5]
    scores = [r["overall_rank"] for r in top5]

    plt.figure(figsize=(8, 5))
    bars = plt.bar(labels, scores, color='skyblue')
    plt.xlabel("Glucose concentration (mM)")
    plt.ylabel("Overall rank (higher = better)")
    plt.title("Top conditions for LPS production (salmonella + clostridium)")
    for bar, score in zip(bars, scores):
        plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                 f"{score:.3f}", ha='center', va='bottom')
    plt.tight_layout()
    plt.savefig(PLOT_OUT, dpi=150)
    print(f"Plot saved to {PLOT_OUT}")
else:
    print("Matplotlib not available – no plot generated")

# ----------------------------------------------------------------------
# Print top conditions to console
# ----------------------------------------------------------------------
print("\nTop 5 conditions (by overall rank):")
for i, r in enumerate(records[:5], 1):
    print(f"{i}. {r['combination']} | {r['ratio']} | Glucose {int(r['glucose_mM'])} mM: "
          f"peak LPS = {r['peak_lps_mM']:.4f} mM at {r['peak_time_h']:.1f} h, "
          f"toxicity_penalty = {r.get('toxicity_penalty', 0.0):.4f}, "
          f"overall rank = {r.get('overall_rank', 0.0):.4f}")
