import cobra
from cobra.io import read_sbml_model
import pandas as pd
from pathlib import Path
from itertools import combinations
from micromol_scorer import create_default_scorer

# ------------------------------------------------------------
# CONFIGURATION
# ------------------------------------------------------------
MODEL_DIR = Path("agora_models")   # change if your models are elsewhere
OUTPUT_CSV = "agora_vaccine_candidates.csv"
NUM_PAIRS = 20                     # test with first 20 pairs only (remove for all)

# ------------------------------------------------------------
# 1. Load all models
# ------------------------------------------------------------
model_files = list(MODEL_DIR.glob("*.sbml")) + list(MODEL_DIR.glob("*.xml"))
print(f"Found {len(model_files)} model files. Loading...")

models = {}
for f in model_files[:10]:   # limit to 10 models for quick test (remove slice for full)
    name = f.stem
    try:
        models[name] = read_sbml_model(f)
        print(f"  Loaded {name}: {len(models[name].reactions)} reactions")
        # Set glucose uptake
        if "EX_glc__D_e" in models[name].reactions:
            models[name].reactions.EX_glc__D_e.lower_bound = -10.0
    except Exception as e:
        print(f"  Failed {name}: {e}")

if len(models) < 2:
    raise ValueError("Need at least 2 models.")

# ------------------------------------------------------------
# 2. Helper
# ------------------------------------------------------------
def get_secreted(model, solution, threshold=1e-6):
    return {r.id for r in model.reactions
            if r.id.startswith("EX_") and solution.fluxes[r.id] > threshold}

# ------------------------------------------------------------
# 3. Pairwise screening
# ------------------------------------------------------------
pairs = list(combinations(models.keys(), 2))
if NUM_PAIRS:
    pairs = pairs[:NUM_PAIRS]
print(f"\nSimulating {len(pairs)} pairs...")

scorer = create_default_scorer()
results = []

for i, (n1, n2) in enumerate(pairs, 1):
    print(f"[{i}/{len(pairs)}] {n1}+{n2}")
    try:
        m1 = models[n1]
        m2 = models[n2]

        sol1 = m1.optimize()
        sol2 = m2.optimize()
        sec1 = get_secreted(m1, sol1)
        sec2 = get_secreted(m2, sol2)

        # Merge models
        comm = m1.copy()
        for rxn in m2.reactions:
            if rxn.id not in comm.reactions:
                comm.add_reaction(rxn.copy())
        bio = [r for r in comm.reactions if "biomass" in r.id.lower()]
        if bio:
            comm.objective = bio[0]
        solc = comm.optimize()
        secc = get_secreted(comm, solc)

        novel = secc - (sec1 | sec2)

        for ex_id in novel:
            rxn = comm.reactions.get_by_id(ex_id)
            met = next((m for m, coeff in rxn.metabolites.items() if abs(coeff) == 1), None)
            met_name = met.name if met and met.name else ex_id
            score = scorer.score_molecule(ex_id, met_name, "", None)
            results.append({
                "Pair": f"{n1}+{n2}",
                "Metabolite": met_name,
                "Vaccine_Score": score.overall_vaccine_score,
                "Recommendation": score.recommendation,
                "Exchange_ID": ex_id
            })
    except Exception as e:
        print(f"  Error: {e}")

# ------------------------------------------------------------
# 4. Save results
# ------------------------------------------------------------
if results:
    df = pd.DataFrame(results).sort_values("Vaccine_Score", ascending=False)
    df.to_csv(OUTPUT_CSV, index=False)
    print(f"\nSaved {len(df)} candidates to {OUTPUT_CSV}")
    print("\nTop 10:")
    print(df.head(10).to_string())
else:
    print("No candidates found.")
    