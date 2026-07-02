#!/usr/bin/env python3
"""
detox_parameter_sweep.py

Sweep `vmax` and `threshold_for_activation` for the detox reaction to find
settings that produce `detoxified_LPS` in the salmonella_detox + clostridium co-culture.
"""
from pathlib import Path
import csv
import json
from micromol_core import (
    create_default_microbes,
    create_metabolite_registry,
    CombinationSimulator,
    create_ecoli_reactions,
    Reaction,
    Microbe,
    Metabolite,
)

OUT_DIR = Path("micromol_results") / "detox_param_sweep"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# parameter ranges
VMAX_VALUES = [5.0, 15.0, 50.0, 150.0]
THRESHOLDS = [0.0, 0.005, 0.01, 0.05]
RATIOS = [(1,1), (1,20)]

# base registries
microbes_base = create_default_microbes()
metabolites = create_metabolite_registry()

# ensure fatty_acids and detoxified_LPS
if "fatty_acids" not in metabolites:
    metabolites["fatty_acids"] = Metabolite(
        id="fatty_acids", name="Fatty acids", smiles="CCCCCCCCCC(=O)O", formula="C10H20O2", mw=172.27, source="nutrient"
    )
if "detoxified_LPS" not in metabolites:
    metabolites["detoxified_LPS"] = Metabolite(
        id="detoxified_LPS", name="Detoxified LPS", smiles="CC(C)C", formula="C8H17O6", mw=234.2, source="metabolite"
    )

results = []

for vmax in VMAX_VALUES:
    for thr in THRESHOLDS:
        # construct detox reaction with current params
        detox_rxn = Reaction(
            id="lipid_A_synthesis",
            name="Detoxified Lipid A synthesis",
            enzyme="LpxL_mutant",
            substrate_ids=["glucose", "phosphate", "fatty_acids"],
            product_ids=["detoxified_LPS"],
            km_values={"glucose": 0.2, "phosphate": 0.1, "fatty_acids": 0.05},
            vmax=float(vmax),
            threshold_for_activation=float(thr)
        )

        # create a fresh microbes dict each iteration
        microbes = create_default_microbes()
        base_rxns = create_ecoli_reactions()
        base_rxns["lipid_A_synthesis"] = detox_rxn
        sal = Microbe(
            id="salmonella_detox",
            name="Salmonella enterica (LPS detox mutant)",
            domain="bacteria",
            reactions=base_rxns,
            native_metabolites={"detoxified_LPS"},
            essential_nutrients={"glucose"},
            max_growth_rate=0.65,
        )
        microbes[sal.id] = sal

        sim = CombinationSimulator(microbes, metabolites)

        for rA, rB in RATIOS:
            biomass_map = {"salmonella_detox": float(rA), "clostridium": float(rB)}
            initial = {"glucose": 100.0, "fatty_acids": 1.0}
            res = sim.simulate_coculture(["salmonella_detox", "clostridium"], initial, duration=48.0, biomass_ratios=biomass_map)
            final = res.get("final_concentrations", {})
            history = res.get("concentration_history", {})
            peak = 0.0
            if history and "detoxified_LPS" in history:
                peak = max(history.get("detoxified_LPS", [0]))

            rec = {
                "vmax": vmax,
                "threshold": thr,
                "ratio": f"{rA}:{rB}",
                "final_detox_mM": final.get("detoxified_LPS", 0.0),
                "peak_detox_mM": peak
            }
            results.append(rec)
            # save per-run
            fname = OUT_DIR / f"run_v{vmax}_t{thr}_r{rA}_{rB}.json"
            with open(fname, "w") as jf:
                json.dump(res, jf, indent=2)

# write summary CSV
csv_path = OUT_DIR / "detox_param_sweep_summary.csv"
with open(csv_path, "w", newline="") as cf:
    writer = csv.DictWriter(cf, fieldnames=["vmax","threshold","ratio","final_detox_mM","peak_detox_mM"]) 
    writer.writeheader()
    for r in results:
        writer.writerow(r)

print(f"Parameter sweep complete. Summary: {csv_path}")
