#!/usr/bin/env python3
"""
debug_detox_monoculture.py

Simulate salmonella_detox alone to debug why lipid_A_synthesis never activates.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from micromol_core import create_default_microbes, create_metabolite_registry, CombinationSimulator
from micromol_core import Microbe, Reaction, Metabolite

# ------------------------------------------------------------
# 1. Build the detox microbe and add missing metabolites
# ------------------------------------------------------------
microbes = create_default_microbes()
metabolites = create_metabolite_registry()

# Ensure fatty_acids is present
if "fatty_acids" not in metabolites:
    metabolites["fatty_acids"] = Metabolite(
        id="fatty_acids",
        name="Fatty acids (pool)",
        smiles="CCCCCCCCCC(=O)O",
        formula="C10H20O2",
        mw=172.27,
        source="nutrient"
    )

# Ensure detoxified_LPS is present
if "detoxified_LPS" not in metabolites:
    metabolites["detoxified_LPS"] = Metabolite(
        id="detoxified_LPS",
        name="Detoxified Lipopolysaccharide",
        smiles="CC(C)CCCC(C)C(=O)OC1CCCC(C1)OC(=O)C2CCCCC2OC(=O)C(C)C",
        formula="C20H36O13",
        mw=468.21,
        source="metabolite"
    )

# Define detox reaction
detox_reaction = Reaction(
    id="lipid_A_synthesis",
    name="Detoxified Lipid A synthesis",
    enzyme="LpxL_mutant",
    substrate_ids=["glucose", "phosphate", "fatty_acids"],
    product_ids=["detoxified_LPS"],
    km_values={"glucose": 0.2, "phosphate": 0.1, "fatty_acids": 0.05},
    vmax=150.0,          # high vmax from sweep
    threshold_for_activation=0.0,   # zero threshold
)

# Create detox microbe (reuse E. coli base reactions)
from micromol_core import create_ecoli_reactions
base_reactions = create_ecoli_reactions()
base_reactions["lipid_A_synthesis"] = detox_reaction

salmonella_detox = Microbe(
    id="salmonella_detox",
    name="Salmonella detox mutant",
    domain="bacteria",
    reactions=base_reactions,
    native_metabolites={"detoxified_LPS"},
    essential_nutrients={"glucose"},
    max_growth_rate=0.65
)

microbes["salmonella_detox"] = salmonella_detox

# ------------------------------------------------------------
# 2. Simulate monoculture
# ------------------------------------------------------------
sim = CombinationSimulator(microbes, metabolites)
initial = {"glucose": 100.0, "fatty_acids": 1.0, "phosphate": 1.0}   # add phosphate explicitly
duration = 48.0

print("Simulating salmonella_detox alone (monoculture)...")
res = sim.simulate_coculture(["salmonella_detox"], initial, duration=duration)

final = res.get("final_concentrations", {})
print("\nFinal concentrations of relevant metabolites:")
for key in ["glucose", "fatty_acids", "phosphate", "detoxified_LPS"]:
    print(f"  {key}: {final.get(key, 0.0):.6f} mM")

# ------------------------------------------------------------
# 3. Check reaction fluxes (if available in result)
# ------------------------------------------------------------
print("\nDEBUGGING NOTE: The current simulator may not expose per-reaction fluxes.")
print("If the detox reaction remains zero, check metabolite IDs and reaction definitions.")
print("You can modify `micromol_core.py` to return reaction fluxes from simulations for deeper inspection.")

print("\nDone.")
