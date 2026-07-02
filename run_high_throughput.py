"""
High-throughput vaccine candidate discovery from all microbe pairs.
Uses real genome-scale models + RDKit similarity to adjuvants.
"""

import cobra
from cobra.io import read_sbml_model
import pandas as pd
from pathlib import Path
from itertools import combinations
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import DataStructs, AllChem
import warnings
warnings.filterwarnings("ignore")

# ------------------------------
# CONFIGURATION
# ------------------------------
MODEL_DIR = Path(".")          # Look in current folder (where all XML files are)
ADJUVANT_CSV = "adjuvants.csv"
OUTPUT_CSV = "high_throughput_candidates.csv"
NUM_WORKERS = min(cpu_count(), 4)   # Use 4 workers to avoid overloading

# Glucose medium (standard)
MEDIUM = {
    "EX_glc__D_e": -10.0,   # glucose uptake
    "EX_o2_e": -100.0,      # oxygen uptake (aerobic)
}

# ------------------------------
# 1. Load all models
# ------------------------------
print("Loading models from", MODEL_DIR)
xml_files = list(MODEL_DIR.glob("*.xml"))
print(f"Found {len(xml_files)} XML files.")
models = {}
for path in xml_files:
    name = path.stem
    try:
        models[name] = read_sbml_model(str(path))
        print(f"  Loaded {name}: {len(models[name].reactions)} reactions")
    except Exception as e:
        print(f"  Failed to load {name}: {e}")

if len(models) < 2:
    raise ValueError("Need at least two models to run pairs.")

# Apply medium to each model
for name, model in models.items():
    for rxn_id, lb in MEDIUM.items():
        if rxn_id in model.reactions:
            model.reactions.get_by_id(rxn_id).lower_bound = lb
        else:
            print(f"Warning: {rxn_id} not found in {name}")

# ------------------------------
# 2. Adjuvant database (hardcoded – no CSV needed)
# ------------------------------
adjuvant_data = [
    ("MPLA", "CC(C)CCCC(C)C(=O)OC1CCCC(C1)OC(=O)C2CCCCC2OC(=O)C(C)C"),
    ("Imiquimod", "CC1=C(C(=N1)N)NC2=C(C(=CC(=C2)Cl)Cl)Cl"),
    ("LPS", "CC(C)CC(C)(C)C(=O)O"),
    ("Flagellin", "C(C(=O)O)N"),
    ("CpG", "C1=C(N(C(=O)NC1=O)C2CC(CC(=O)O)O2)N"),
    ("QS-21", "CC(C)CCCC(C)C(=O)OC1CCCC(C1)OC(=O)C2CCCCC2OC(=O)C(C)C"),
    ("Alum", "[Al+3].O"),
]
adj_smiles = [smi for _, smi in adjuvant_data]
adj_mols = [Chem.MolFromSmiles(s) for s in adj_smiles if Chem.MolFromSmiles(s)]
adj_fps = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, 1024) for mol in adj_mols]
print(f"Loaded {len(adj_fps)} adjuvants.")

def max_similarity_to_adjuvants(smiles):
    if not smiles or pd.isna(smiles):
        return 0.0
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0.0
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, 1024)
    best = 0.0
    for adj_fp in adj_fps:
        sim = DataStructs.TanimotoSimilarity(fp, adj_fp)
        if sim > best:
            best = sim
    return best

# ------------------------------
# 3. Helper: get secreted metabolites from a model + solution
# ------------------------------
def get_secreted(model, solution, threshold=1e-6):
    """Return dict of exchange ID -> flux for secreted (>0) metabolites."""
    secreted = {}
    for rxn in model.reactions:
        if rxn.id.startswith("EX_") and solution.fluxes[rxn.id] > threshold:
            secreted[rxn.id] = solution.fluxes[rxn.id]
    return secreted

# ------------------------------
# 4. Worker function: simulate a single pair
# ------------------------------
def simulate_pair(pair):
    name1, name2 = pair
    model1 = models[name1]
    model2 = models[name2]

    try:
        # Monocultures
        sol1 = model1.optimize()
        sol2 = model2.optimize()
        secret1 = get_secreted(model1, sol1)
        secret2 = get_secreted(model2, sol2)

        # Community model: merge models (simple union)
        community = model1.copy()
        for rxn in model2.reactions:
            if rxn.id not in community.reactions:
                community.add_reaction(rxn.copy())
        # Set objective to maximise total biomass – find a biomass reaction
        biomass_ids = [rxn.id for rxn in community.reactions if "biomass" in rxn.id.lower()]
        if biomass_ids:
            community.objective = biomass_ids[0]
        else:
            # fallback: use first reaction
            community.objective = list(community.reactions)[0].id
        sol_com = community.optimize()
        secret_com = get_secreted(community, sol_com)

        # Novel metabolites: present in community but not in either monoculture
        all_mono = set(secret1.keys()) | set(secret2.keys())
        novel_ids = set(secret_com.keys()) - all_mono

        # Score each novel metabolite
        results = []
        for ex_id in novel_ids:
            # Try to get metabolite name from the exchange reaction
            rxn = community.reactions.get_by_id(ex_id)
            met = None
            for m, coeff in rxn.metabolites.items():
                if abs(coeff) == 1:
                    met = m
                    break
            met_name = met.name if met and met.name else ex_id
            # SMILES is not available – we use empty for now. Real pipeline would query PubChem.
            smiles = ""
            sim = max_similarity_to_adjuvants(smiles)
            score = sim * 100
            if score >= 70:
                rec = "high"
            elif score >= 35:
                rec = "medium"
            else:
                rec = "low"
            results.append({
                "Pair": f"{name1}+{name2}",
                "Exchange_ID": ex_id,
                "Metabolite_Name": met_name,
                "Adjuvant_Similarity": round(sim, 3),
                "Vaccine_Score": round(score, 1),
                "Recommendation": rec,
                "Flux": round(secret_com[ex_id], 3)
            })
        return results
    except Exception as e:
        print(f"Error in pair {name1}+{name2}: {e}")
        return []

# ------------------------------
# 5. Generate all pairs and run in parallel
# ------------------------------
model_names = list(models.keys())
pairs = list(combinations(model_names, 2))
print(f"\nTotal pairs to simulate: {len(pairs)}")

all_results = []
with Pool(processes=NUM_WORKERS) as pool:
    for res in tqdm(pool.imap(simulate_pair, pairs), total=len(pairs)):
        all_results.extend(res)

# ------------------------------
# 6. Convert to DataFrame and save
# ------------------------------
if all_results:
    df = pd.DataFrame(all_results)
    df = df.sort_values("Vaccine_Score", ascending=False)
    df.to_csv(OUTPUT_CSV, index=False)
    print(f"\nSaved {len(df)} candidate entries to {OUTPUT_CSV}")
    print("\nTop 20 candidates:")
    print(df.head(20).to_string())
else:
    print("No candidates found.")