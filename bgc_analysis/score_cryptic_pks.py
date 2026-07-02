#!/usr/bin/env python3
"""
score_cryptic_pks.py

Generate plausible polyketide scaffolds for cryptic T1PKS cluster (Region 1.6)
and score them for vaccine potential.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from micromol_scorer import create_default_scorer

scorer = create_default_scorer()

candidates = [
    {
        "id": "pks_macrolactone_14",
        "name": "14-membered macrolactone",
        "smiles": "CC1CCCC(=O)OC(C)CC(C)CC(C)CC1=O",
        "mw": 300
    },
    {
        "id": "pks_macrolactone_16",
        "name": "16-membered macrolactone",
        "smiles": "CC1CCCCCC(=O)OC(C)CC(C)CC(C)CC1=O",
        "mw": 340
    },
    {
        "id": "pks_linear_polyene",
        "name": "Linear polyene",
        "smiles": "CC=CC=CC=CC=CC(=O)O",
        "mw": 220
    },
    {
        "id": "pks_aromatic",
        "name": "Aromatic polyketide (tetracene)",
        "smiles": "CC1=C2C(=C3C(=C1O)C(=O)C4=CC=CC=C4C3=O)O",
        "mw": 320
    },
    {
        "id": "pks_ansamycin_like",
        "name": "Ansamycin-like macrolactam",
        "smiles": "CC1CCCC(=O)NCC(C)CC(C)CC1=O",
        "mw": 310
    },
    {
        "id": "pks_spiroketal",
        "name": "Spiroketal polyketide",
        "smiles": "CC1CC2(C(CC1=O)OC(O2)C)C",
        "mw": 280
    },
    {
        "id": "pks_linear_polyketide",
        "name": "Linear polyketide with beta-hydroxy",
        "smiles": "CCCC(=O)CC(O)C=CCCC(=O)O",
        "mw": 260
    },
    {
        "id": "pks_terpenoid_hybrid",
        "name": "Terpenoid-polyketide hybrid",
        "smiles": "CC1=CCC(C(C)C)CC1C(=O)O",
        "mw": 250
    }
]

print("Scoring cryptic PKS scaffolds from Region 1.6 (Streptomyces avermitilis)\n")
results = []
for cand in candidates:
    score = scorer.score_molecule(
        molecule_id=cand["id"],
        molecule_name=cand["name"],
        smiles=cand["smiles"],
        molecular_weight=cand["mw"]
    )
    results.append({
        "id": cand["id"],
        "Name": cand["name"],
        "Vaccine_Score": getattr(score, 'overall_vaccine_score', None),
        "Recommendation": getattr(score, 'recommendation', ''),
        "Toxicity": getattr(score, 'toxicity_risk', None),
        "Immunogenicity": getattr(score, 'immunogenicity_score', None)
    })

results_sorted = sorted(results, key=lambda x: (x['Vaccine_Score'] or 0), reverse=True)

print("Ranked candidates:\n")
for i, r in enumerate(results_sorted, 1):
    print(f"{i}. {r['Name']}: Score = {r['Vaccine_Score']:.1f} ({r['Recommendation']}), Toxicity = {r['Toxicity']:.3f}")

best = results_sorted[0]
print("\n" + "="*50)
print(f"TOP CANDIDATE from cryptic T1PKS (Region 1.6): {best['Name']}")
print(f"Vaccine score: {best['Vaccine_Score']:.1f} / 100 – {best['Recommendation']}")
print("="*50)
