#!/usr/bin/env python3
"""
produce_publication_summary.py

Post-process optimisation runs to compute integrated LPS yield
and recommend harvest times. Produces:
 - micromol_results/optimisation/glucose_biomass_sweep/harvest_summary.csv
 - micromol_results/optimisation/glucose_biomass_sweep/harvest_heatmap.png

Usage:
  python bgc_analysis/produce_publication_summary.py
"""
from pathlib import Path
import csv
import json
import numpy as np
import math

ROOT = Path("micromol_results") / "optimisation" / "glucose_biomass_sweep"
IN_CSV = ROOT / "optimisation_summary.csv"
OUT_CSV = ROOT / "harvest_summary.csv"
HEATMAP_PNG = ROOT / "harvest_heatmap.png"

TARGET = "polysaccharide_lps"
TOXIC = ["ethanol", "lactate", "acetate"]

try:
    import matplotlib.pyplot as plt
    HAS_PLT = True
except Exception:
    HAS_PLT = False

if not IN_CSV.exists():
    raise SystemExit(f"Input summary not found: {IN_CSV}")

rows = []
with open(IN_CSV, "r", newline="") as f:
    reader = csv.DictReader(f)
    for r in reader:
        rows.append(r)

out_rows = []
for r in rows:
    combo = r["combination"]
    glucose = int(float(r["glucose_mM"]))
    ratio = r["ratio"]
    # reconstruct run dir used earlier: combo_safe and ratio string with underscore
    combo_safe = combo.replace(" ", "_").replace("+", "_")
    ratio_undersc = ratio.replace(":", "_")
    run_dir = Path("micromol_results") / "optimisation" / combo_safe / f"glucose_{glucose}_ratio_{ratio_undersc}"
    jd = run_dir / "micromol_detailed.json"
    if not jd.exists():
        # skip missing
        continue
    data = json.loads(jd.read_text())
    if isinstance(data, list) and data:
        data = data[0]

    time = data.get("time_steps", [])
    hist = data.get("concentration_history", {})
    lps = hist.get(TARGET, [])
    # compute integrated yield over first 2 hours (or full if shorter)
    if not time or not lps:
        integrated = 0.0
        peak = 0.0
        harvest_time = 0.0
    else:
        max_time = 2.0
        # select indices up to max_time
        idx = [i for i,t in enumerate(time) if t <= max_time]
        if not idx:
            # fallback: use full
            xs = np.array(time)
            ys = np.array(lps)
        else:
            xs = np.array([time[i] for i in idx])
            ys = np.array([lps[i] for i in idx])

        # trapezoidal integration (manual to avoid numpy.trapz absence)
        if xs.size < 2:
            integrated = float(0.0)
        else:
            dx = xs[1:] - xs[:-1]
            mid = (ys[1:] + ys[:-1]) / 2.0
            integrated = float((mid * dx).sum())
        peak = float(np.max(ys))
        # recommended harvest: first time reaching 95% of peak (within window)
        thresh = 0.95 * peak if peak > 0 else 0.0
        harvest_time = 0.0
        if thresh > 0:
            for tval, yval in zip(xs, ys):
                if yval >= thresh:
                    harvest_time = float(tval)
                    break

    # toxicity: use max of toxic metabolites over same window
    tox_vals = []
    for tox in TOXIC:
        hist_t = hist.get(tox, [])
        if hist_t and xs.size>0:
            # align lengths: hist arrays are same length as full time, so pick indices
            tox_slice = [hist_t[t_idx] for t_idx in range(min(len(hist_t), len(data.get("time_steps", []))))]
            tox_vals.append(float(np.max(tox_slice)))
        elif hist_t:
            tox_vals.append(float(np.max(hist_t)))
        else:
            tox_vals.append(0.0)

    toxicity_score = float(max(tox_vals) if tox_vals else 0.0)

    out_rows.append({
        "combination": combo,
        "microbes": r.get("microbes", ""),
        "glucose_mM": glucose,
        "ratio": ratio,
        "integrated_lps_mM_h": integrated,
        "peak_lps_mM": peak,
        "harvest_time_h": harvest_time,
        "toxicity_max_mM": toxicity_score
    })

# write CSV
OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
with open(OUT_CSV, "w", newline="") as f:
    fieldnames = ["combination","microbes","glucose_mM","ratio","integrated_lps_mM_h","peak_lps_mM","harvest_time_h","toxicity_max_mM"]
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    for rr in out_rows:
        writer.writerow(rr)

print(f"Harvest summary saved to {OUT_CSV}")

if HAS_PLT and out_rows:
    # pivot for heatmap: x=glucose, y=ratio, value=integrated
    ratios = sorted(list({rr["ratio"] for rr in out_rows}), key=lambda s: (int(s.split(":")[0]), int(s.split(":")[1])))
    glucoses = sorted(list({int(rr["glucose_mM"]) for rr in out_rows}))
    grid = np.zeros((len(ratios), len(glucoses)))
    for rr in out_rows:
        i = ratios.index(rr["ratio"])
        j = glucoses.index(int(rr["glucose_mM"]))
        grid[i, j] = rr["integrated_lps_mM_h"]

    import matplotlib.pyplot as plt
    plt.figure(figsize=(6,4))
    im = plt.imshow(grid, aspect='auto', cmap='viridis')
    plt.colorbar(im, label='Integrated LPS (mM·h)')
    plt.yticks(range(len(ratios)), ratios)
    plt.xticks(range(len(glucoses)), [str(g) for g in glucoses])
    plt.xlabel('Glucose (mM)')
    plt.ylabel('Ratio (A:B)')
    plt.title('Integrated LPS (0-2 h)')
    plt.tight_layout()
    plt.savefig(HEATMAP_PNG, dpi=200)
    print(f"Heatmap saved to {HEATMAP_PNG}")
