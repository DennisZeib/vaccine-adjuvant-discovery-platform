#!/usr/bin/env python3
"""
test_detox_with_fattyacids.py

Add fatty acids to the initial medium for salmonella_detox + clostridium,
then simulate and check for detoxified_LPS production.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

from micromol_core import create_default_microbes, create_metabolite_registry, CombinationSimulator

# Add the detox microbe and new metabolite to the existing registries
from micromol_core import Microbe, Reaction, Metabolite

# Get existing registries
microbes = create_default_microbes()
metabolites = create_metabolite_registry()

# Add fatty_acids if not already present (it is needed)
if "fatty_acids" not in metabolites:
    metabolites["fatty_acids"] = Metabolite(
        id="fatty_acids",
        name="Fatty acids (pool)",
        smiles="CCCCCCCCCC(=O)O",
        formula="C10H20O2",
        mw=172.27,
        source="nutrient"
    )

# Add detoxified_LPS as a new metabolite
if "detoxified_LPS" not in metabolites:
    metabolites["detoxified_LPS"] = Metabolite(
        id="detoxified_LPS",
        name="Detoxified Lipopolysaccharide",
        smiles="CC(C)CCCC(C)C(=O)OC1CCCC(C1)OC(=O)C2CCCCC2OC(=O)C(C)C",
        formula="C20H36O13",
        mw=468.21,
        source="metabolite"
    )

# Define the detox microbe reaction
detox_reaction = Reaction(
    id="lipid_A_synthesis",
    name="Detoxified Lipid A synthesis",
    enzyme="LpxL_mutant",
    substrate_ids=["glucose", "phosphate", "fatty_acids"],
    product_ids=["detoxified_LPS"],
    km_values={"glucose": 0.2, "phosphate": 0.1, "fatty_acids": 0.05},
    vmax=15.0,
    threshold_for_activation=0.01
)

# Reuse E. coli reactions
from micromol_core import create_ecoli_reactions
base_reactions = create_ecoli_reactions()
base_reactions["lipid_A_synthesis"] = detox_reaction

salmonella_detox = Microbe(
    id="salmonella_detox",
    name="Salmonella enterica (LPS detox mutant)",
    domain="bacteria",
    reactions=base_reactions,
    native_metabolites={"detoxified_LPS"},
    essential_nutrients={"glucose"},
    max_growth_rate=0.65
)

# Replace or add to microbes dict
microbes["salmonella_detox"] = salmonella_detox

# Now simulate co-culture with clostridium
sim = CombinationSimulator(microbes, metabolites)

# Initial nutrients: glucose + fatty acids
initial = {
    "glucose": 100.0,
    "fatty_acids": 1.0
}

print("Simulating salmonella_detox + clostridium with fatty acids added...")
res = sim.simulate_coculture(["salmonella_detox", "clostridium"], initial, duration=48.0)

final = res.get("final_concentrations", {})
detox_yield = final.get("detoxified_LPS", 0.0)
print(f"\nDetoxified LPS final concentration: {detox_yield:.6f} mM")

if detox_yield > 0:
    print("SUCCESS: detoxified_LPS produced!")
else:
    print("Still zero. Now we need a parameter sweep (Option C).")
