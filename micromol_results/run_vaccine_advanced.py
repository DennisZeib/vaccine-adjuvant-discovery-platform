"""
Advanced vaccine discovery pipeline using AGORA2 microbial models.
Simulates pairwise co‑cultures, detects novel secreted metabolites,
and scores them for vaccine potential.
"""

import cobra
from cobra.io import read_sbml_model
import pandas as pd
from pathlib import Path
from itertools import combinations
from micromol_scorer import create_default_scorer
import time
import sys

# Increase recursion limit to avoid "maximum recursion depth exceeded" errors
sys.setrecursionlimit(10000)

# ============================================================
# CONFIGURATION – adjust these as needed
# ============================================================
MODEL_DIR = Path("agora_models")          # folder containing all .xml models
OUTPUT_CSV = "vaccine_candidates_advanced.csv"
NUM_MODELS = 20                           # start with 20 models (190 pairs)
# Set NUM_MODELS = None to use all models in the folder (not recommended for >100 models)

# ============================================================
# 1. Load models
# ============================================================
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
        # Set glucose uptake if exchange exists
        if "EX_glc__D_e" in models[name].reactions:
            models[name].reactions.EX_glc__D_e.lower_bound = -10.0
    except Exception as e:
        print(f"FAILED: {e}")

if len(models) < 2:
    raise ValueError("Need at least 2 successfully loaded models.")

# ============================================================
# 2. Helper function to get secreted metabolites
# ============================================================
def get_secreted(model, solution, threshold=1e-6):
    """Return set of exchange reaction IDs with positive flux."""
    return {rxn.id for rxn in model.reactions
            if rxn.id.startswith("EX_") and solution.fluxes[rxn.id] > threshold}

# ============================================================
# 3. Pairwise simulation
# ============================================================
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

        # Monoculture secretions
        sol1 = m1.optimize()
        sol2 = m2.optimize()
        sec1 = get_secreted(m1, sol1)
        sec2 = get_secreted(m2, sol2)

        # Build community model by merging
        community = m1.copy()
        for rxn in m2.reactions:
            if rxn.id not in community.reactions:
                # Use add_reactions (plural) with a list – works in all COBRApy versions
                community.add_reactions([rxn.copy()])

        # Set objective to a biomass reaction if present
        biomass_rxns = [r for r in community.reactions if "biomass" in r.id.lower()]
        if biomass_rxns:
            community.objective = biomass_rxns[0]

        # Optimise community
        solc = community.optimize()
        secc = get_secreted(community, solc)

        # Novel = produced only in co‑culture
        novel = secc - (sec1 | sec2)

        # Score each novel metabolite
        for ex_id in novel:
            try:
                rxn = community.reactions.get_by_id(ex_id)
                # Find the exchanged metabolite (coefficient ±1)
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
                # Skip single metabolite if it causes error
                print(f"    Warning: could not score {ex_id}: {e}")

        elapsed = time.time() - start_time
        print(f"    Done in {elapsed:.1f}s, found {len(novel)} novel metabolites")

    except RecursionError as e:
        print(f"    RecursionError: {e} – skipping this pair")
    except Exception as e:
        print(f"    Error: {e} – skipping")

# ============================================================
# 4. Save and display results
# ============================================================
if results:
    df = pd.DataFrame(results)
    df = df.sort_values("Vaccine_Score", ascending=False)
    df.to_csv(OUTPUT_CSV, index=False)
    print(f"\n✅ Saved {len(df)} candidate entries to {OUTPUT_CSV}")
    print("\n🏆 TOP 20 VACCINE CANDIDATES:")
    print(df.head(20).to_string(index=False))
else:
    print("\n❌ No candidates found. Try running with more models or different medium conditions.")

print("\nDone.")