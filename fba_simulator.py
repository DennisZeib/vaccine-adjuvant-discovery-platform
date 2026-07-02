from cobra.io import read_sbml_model
from cobra.flux_analysis import pfba
import pandas as pd

# Load real models
ecoli = read_sbml_model("models/iML1515.xml")
yeast = read_sbml_model("models/iND750.xml")

# Set glucose as the only carbon source
ecoli.reactions.EX_glc__D_e.lower_bound = -10.0  # uptake 10 mmol/gDW/h
yeast.reactions.EX_glc__D_e.lower_bound = -10.0

# Simulate each alone
sol_ecoli = pfba(ecoli)
sol_yeast = pfba(yeast)

# Extract produced metabolites (exchange reactions with positive flux)
ecoli_produced = [r.id for r in ecoli.reactions if r.id.startswith("EX_") and sol_ecoli.fluxes[r.id] > 0]
yeast_produced = [r.id for r in yeast.reactions if r.id.startswith("EX_") and sol_yeast.fluxes[r.id] > 0]

# Co‑culture: simple merge (but real community FBA is more complex)
# For now, just take the union of produced metabolites
co_produced = set(ecoli_produced) | set(yeast_produced)

# Novel molecules = produced in co‑culture but not in either alone
novel = co_produced - (set(ecoli_produced) & set(yeast_produced))
print("Novel metabolites in co‑culture:", novel)