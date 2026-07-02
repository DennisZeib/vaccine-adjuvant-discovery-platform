#!/usr/bin/env python3
"""
detox_glucose_sweep.py

Run a glucose concentration sweep for `salmonella_detox` monoculture and record
`detoxified_LPS` peak, final and integrated (0-2 h) yields.
"""
from pathlib import Path
import csv
import json
import math
from micromol_core import (
    create_default_microbes,
    create_metabolite_registry,
    CombinationSimulator,
    Microbe,
    Reaction,
    Metabolite,
    create_ecoli_reactions,
)

OUT_DIR = Path("micromol_results") / "optimisation" / "detox_glucose_sweep"
OUT_DIR.mkdir(parents=True, exist_ok=True)

GLU_VALUES = [10.0, 20.0, 50.0, 100.0]
duration = 48.0

# prepare microbe + metabolites
microbes = create_default_microbes()
metabolites = create_metabolite_registry()

if "fatty_acids" not in metabolites:
    metabolites["fatty_acids"] = Metabolite(
        id="fatty_acids", name="Fatty acids", smiles="CCCCCCCCCC(=O)O", formula="C10H20O2", mw=172.27, source="nutrient"
    )
if "detoxified_LPS" not in metabolites:
    metabolites["detoxified_LPS"] = Metabolite(
        id="detoxified_LPS", name="Detoxified LPS", smiles="CC(C)C", formula="C8H17O6", mw=234.2, source="metabolite"
    )

detox_rxn = Reaction(
    id="lipid_A_synthesis",
    name="Detoxified Lipid A synthesis",
    enzyme="LpxL_mutant",
    substrate_ids=["glucose", "phosphate", "fatty_acids"],
    product_ids=["detoxified_LPS"],
    km_values={"glucose": 0.2, "phosphate": 0.1, "fatty_acids": 0.05},
    vmax=150.0,
    threshold_for_activation=0.0,
)

base_rxns = create_ecoli_reactions()
base_rxns["lipid_A_synthesis"] = detox_rxn

sal = Microbe(
    id="salmonella_detox",
    name="Salmonella detox mutant",
    domain="bacteria",
    reactions=base_rxns,
    native_metabolites={"detoxified_LPS"},
    essential_nutrients={"glucose"},
    max_growth_rate=0.65,
)

microbes[sal.id] = sal

sim = CombinationSimulator(microbes, metabolites)

rows = []
for g in GLU_VALUES:
    initial = {"glucose": float(g), "fatty_acids": 1.0, "phosphate": 1.0}
    res = sim.simulate_coculture(["salmonella_detox"], initial, duration=duration)
    final = res.get("final_concentrations", {})
    history = res.get("concentration_history", {})
    times = res.get("time_points", [])

    detox_hist = history.get("detoxified_LPS", []) if history else []
    peak = max(detox_hist) if detox_hist else 0.0
    final_val = final.get("detoxified_LPS", 0.0)

    # integrated 0-2h using simple trapezoid if time_points present
    integ = 0.0
    if times and detox_hist:
        # collect indices up to 2.0
        pairs = [(t, c) for t, c in zip(times, detox_hist) if t <= 2.0]
        if len(pairs) >= 2:
            xs = [p[0] for p in pairs]
            ys = [p[1] for p in pairs]
            for i in range(1, len(xs)):
                integ += 0.5 * (ys[i] + ys[i-1]) * (xs[i] - xs[i-1])
    rows.append({"glucose_mM": g, "peak_mM": peak, "final_mM": final_val, "integrated_0_2h_mM_h": integ})

    # save per-run
    with open(OUT_DIR / f"detox_glucose_{int(g)}_detailed.json", "w") as jf:
        json.dump(res, jf, indent=2)

# write summary CSV
csv_path = OUT_DIR / "detox_glucose_sweep_summary.csv"
with open(csv_path, "w", newline="") as cf:
    writer = csv.DictWriter(cf, fieldnames=["glucose_mM","peak_mM","final_mM","integrated_0_2h_mM_h"]) 
    writer.writeheader()
    for r in rows:
        writer.writerow(r)

print(f"Detox glucose sweep complete. Summary: {csv_path}")
