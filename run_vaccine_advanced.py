import cobra
from cobra.io import read_sbml_model
import pandas as pd
from pathlib import Path
from itertools import combinations
from micromol_scorer import create_default_scorer
import time
import sys

sys.setrecursionlimit(10000)

MODEL_DIR = Path("agora_models")
OUTPUT_CSV = "vaccine_candidates_advanced.csv"
NUM_MODELS = none

print("Scanning for model files...")
model_files = list(MODEL_DIR.glob("*.xml"))
if not model_files:
    raise FileNotFoundError(f"No .xml files found in {MODEL_DIR}")

if NUM_MODELS and NUM_MODELS > 0:
    model_files = model_files[:NUM_MODELS]
print(f"Will load {len(model_files)} models.")

models = {}
for f in model_files:
    name = f.stem
    try:
        print(f"Loading {name} ...", end=" ", flush=True)
        models[name] = read_sbml_model(f)
        print(f"OK ({len(models[name].reactions)} reactions)")
        if "EX_glc__D_e" in models[name].reactions:
            models[name].reactions.EX_glc__D_e.lower_bound = -10.0
    except Exception as e:
        print(f"FAILED: {e}")

if len(models) < 2:
    raise ValueError("Need at least 2 successfully loaded models.")

def get_secreted(model, solution, threshold=1e-6):
    return {rxn.id for rxn in model.reactions
            if rxn.id.startswith("EX_") and solution.fluxes[rxn.id] > threshold}

pairs = list(combinations(models.keys(), 2))
print(f"\nTotal pairs to simulate: {len(pairs)}")
scorer = create_default_scorer()
results = []

for i, (n1, n2) in enumerate(pairs, 1):
    start_time = time.time()
    print(f"[{i}/{len(pairs)}] {n1} + {n2}")
    try:
        m1 = models[n1]
        m2 = models[n2]
        sol1 = m1.optimize()
        sol2 = m2.optimize()
        sec1 = get_secreted(m1, sol1)
        sec2 = get_secreted(m2, sol2)
        community = m1.copy()
        for rxn in m2.reactions:
            if rxn.id not in community.reactions:
                community.add_reactions([rxn.copy()])
        biomass_rxns = [r for r in community.reactions if "biomass" in r.id.lower()]
        if biomass_rxns:
            community.objective = biomass_rxns[0]
        solc = community.optimize()
        secc = get_secreted(community, solc)
        novel = secc - (sec1 | sec2)
        for ex_id in novel:
            try:
                rxn = community.reactions.get_by_id(ex_id)
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
                print(f"    Warning: could not score {ex_id}: {e}")
        elapsed = time.time() - start_time
        print(f"    Done in {elapsed:.1f}s, found {len(novel)} novel metabolites")
    except RecursionError as e:
        print(f"    RecursionError: {e} – skipping this pair")
    except Exception as e:
        print(f"    Error: {e} – skipping")

if results:
    df = pd.DataFrame(results).sort_values("Vaccine_Score", ascending=False)
    df.to_csv(OUTPUT_CSV, index=False)
    print(f"\nSaved {len(df)} candidate entries to {OUTPUT_CSV}")
    print("\nTOP 20 VACCINE CANDIDATES:")
    print(df.head(20).to_string(index=False))
else:
    print("\nNo candidates found.")
