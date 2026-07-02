#!/usr/bin/env python3
"""
compare_toxicity.py

Compare toxicity and immunogenicity of wild-type LPS vs detoxified_LPS.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from micromol_scorer import create_default_scorer

scorer = create_default_scorer()

lps_score = scorer.score_molecule(
    molecule_id="polysaccharide_lps",
    molecule_name="Lipopolysaccharide (wild-type)",
    smiles="CC(C)CC(C)(C)C(=O)O",
    molecular_weight=468.21
)

detox_score = scorer.score_molecule(
    molecule_id="detoxified_LPS",
    molecule_name="Detoxified Lipopolysaccharide",
    smiles="CC(C)CCCC(C)C(=O)OC1CCCC(C1)OC(=O)C2CCCCC2OC(=O)C(C)C",
    molecular_weight=468.21
)

print("="*50)
print("TOXICITY & IMMUNOGENICITY COMPARISON")
print("="*50)
print(f"Wild-type LPS   : Toxicity risk = {lps_score.toxicity_risk:.3f}, Immunogenicity = {lps_score.immunogenicity_score:.3f}")
print(f"Detoxified LPS  : Toxicity risk = {detox_score.toxicity_risk:.3f}, Immunogenicity = {detox_score.immunogenicity_score:.3f}")
print("="*50)

if detox_score.toxicity_risk < lps_score.toxicity_risk:
    print("✅ Detoxified LPS is predicted to be LESS TOXIC than wild-type LPS.")
else:
    print("⚠️ Detoxified LPS toxicity not lower – adjust reaction or structure.")

if detox_score.immunogenicity_score > lps_score.immunogenicity_score:
    print("✅ Detoxified LPS is predicted to be MORE IMMUNOGENIC.")
else:
    print("ℹ️ Detoxified LPS has similar or lower immunogenicity – may still be useful as a safer adjuvant.")
