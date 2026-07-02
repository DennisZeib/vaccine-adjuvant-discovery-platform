#!/usr/bin/env python3
"""
score_filipin.py

Score filipin (polyene macrolide) with the micromol scorer.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from micromol_scorer import create_default_scorer

scorer = create_default_scorer()

# Filipin SMILES (simplified placeholder)
smiles = "CC(C)CC1=CC(=O)OC2CC(CC(C2C(=O)C=C(C(=O)O1)C)O)C"
score = scorer.score_molecule("filipin", "Filipin (polyene macrolide)", smiles, molecular_weight=700)

print(f"Filipin vaccine score: {score.overall_vaccine_score:.1f} ({score.recommendation})")
print(f"Toxicity risk: {score.toxicity_risk:.3f}")
