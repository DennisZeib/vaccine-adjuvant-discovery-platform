#!/usr/bin/env python3
"""
Run a detox mutant simulation: defines a `salmonella_detox` microbe and
simulates with `clostridium` to check detoxified_LPS production.
"""
from pathlib import Path
import json

from micromol_core import (
    create_default_microbes,
    create_metabolite_registry,
    CombinationSimulator,
    Reaction,
    Microbe,
    create_ecoli_reactions,
    Metabolite,
)

OUT_DIR = Path("micromol_results") / "detox_runs"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# setup
microbes = create_default_microbes()
metabolites = create_metabolite_registry()
sim = CombinationSimulator(microbes, metabolites)

# add new metabolite: detoxified_LPS
if "detoxified_LPS" not in metabolites:
    metabolites["detoxified_LPS"] = Metabolite(
        id="detoxified_LPS",
        name="Detoxified LPS (virtual)",
        smiles="CC(=O)OCCOP(=O)(O)O",  # placeholder
        formula="C8H17O6P",
        mw=234.2,
        source="metabolite",
    )

# base reactions from E. coli to reuse simple reactions
base_rxns = create_ecoli_reactions()

# define lipid_A_synthesis (detoxified)
lipid_rxn = Reaction(
    id="lipid_A_synthesis",
    name="Detoxified Lipid A synthesis",
    enzyme="LpxL_mutant",
    substrate_ids=["glucose", "ATP"],
    product_ids=["detoxified_LPS"],
    km_values={"glucose": 0.2, "ATP": 0.1},
    vmax=15.0,
)

flag_rxn = Reaction(
    id="flagellin_synthesis",
    name="Flagellin synthesis",
    enzyme="Flagellin synthase",
    substrate_ids=["alanine", "ATP"],
    product_ids=["flagellin"],
    km_values={"alanine": 0.5, "ATP": 0.1},
    vmax=8.0,
)

# create reaction dict for detox microbe
detox_reactions = dict(base_rxns)
detox_reactions.update({"lipid_A_synthesis": lipid_rxn, "flagellin_synthesis": flag_rxn})

# define microbe
sal_detox = Microbe(
    id="salmonella_detox",
    name="Salmonella enterica (LPS detox mutant)",
    domain="bacteria",
    reactions=detox_reactions,
    native_metabolites={"detoxified_LPS", "flagellin"},
    essential_nutrients={"glucose"},
    max_growth_rate=0.65,
)

# add to microbes dict
microbes[sal_detox.id] = sal_detox

# ensure sim knows about new microbe and metabolite registry
sim = CombinationSimulator(microbes, metabolites)

# run small sweep: glucose 100 mM, ratios 1:1 and 1:20 with clostridium
glucose = 100.0
ratios = [(1, 1), (1, 20)]
results = []
for rA, rB in ratios:
    biomass_map = {"salmonella_detox": float(rA), "clostridium": float(rB)}
    initial = {"glucose": float(glucose)}
    print(f"Running detox mutant | glucose={glucose} mM | ratio={rA}:{rB} ...")
    res = sim.simulate_coculture(["salmonella_detox", "clostridium"], initial, duration=48.0, biomass_ratios=biomass_map)
    out = {
        "ratio": f"{rA}:{rB}",
        "final_concentrations": res.get("final_concentrations", {}),
        "peak_lps": max(res.get("concentration_history", {}).get("detoxified_LPS", [0])) if res.get("concentration_history") else 0.0,
        "time_steps": res.get("time_steps", []),
    }
    results.append(out)
    # save detailed
    with open(OUT_DIR / f"detox_glucose_{int(glucose)}_ratio_{rA}_{rB}.json", "w") as jf:
        json.dump(res, jf, indent=2)

with open(OUT_DIR / "detox_summary.json", "w") as jf:
    json.dump(results, jf, indent=2)

print("Detox runs complete. Summary:")
for r in results:
    fc = r["final_concentrations"]
    print(f"Ratio {r['ratio']}: detoxified_LPS final = {fc.get('detoxified_LPS', 0.0):.4f} mM, peak = {r['peak_lps']:.4f}")
