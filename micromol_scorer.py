"""
MicroMol Vaccine Potential Scorer
Evaluates molecules for vaccine/adjuvant potential using manual heuristics (no RDKit required).
"""

import logging
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass

# RDKit is optional; we don't use it (manual scoring only)
try:
    from rdkit import Chem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

logger = logging.getLogger(__name__)


@dataclass
class VaccinePotentialScore:
    """Represents vaccine potential assessment of a molecule."""
    molecule_id: str
    molecule_name: str
    immunogenicity_score: float  # 0-1
    novelty_score: float         # 0-1
    toxicity_risk: float         # 0-1 (lower better)
    similarity_to_adjuvants: float
    similarity_to_toxins: float
    overall_vaccine_score: float  # 0-100
    recommendation: str           # "high", "medium", "low", "excluded"
    rationale: List[str]


class KnownMoleculeDatabase:
    """Built‑in reference of adjuvants and toxins (for completeness)."""
    VACCINE_ADJUVANTS = {
        "mpla": {"smiles": "CC(C)CCCC(C)C(=O)OC1CCCC(C1)OC(=O)C2CCCCC2OC(=O)C(C)C", "name": "MPLA"},
        "lps": {"smiles": "CC(C)CC(C)(C)C(=O)O", "name": "LPS"},
        "flagellin": {"smiles": "C(C(=O)O)N", "name": "Flagellin"},
        "cpG": {"smiles": "C1=C(N(C(=O)NC1=O)C2CC(CC(=O)O)O2)N", "name": "CpG"},
    }
    KNOWN_TOXINS = {
        "aflatoxin_b1": {"smiles": "CC1=C2C(=C3C(=C1OC)C(=O)OC3=CC4=CO4)OC2=O", "name": "Aflatoxin B1"},
    }


class VaccineScorer:
    """Scoring engine that works without RDKit using manual molecule database."""

    def __init__(self):
        self.db = KnownMoleculeDatabase()
        # Manual similarity scores for molecules the simulator produces
        self.manual_adjuvant_sim = {
            "polysaccharide_lps": 0.85,
            "flagellin": 0.90,
            "mannan_polysaccharide": 0.75,
            "beta_glucan": 0.80,
            "viral_capsid": 0.65,
            "viral_particle_with_spike": 0.70,
            "outer_membrane_vesicle": 0.88,
            "ATP": 0.05,
            "glycerol": 0.02,
            "ethanol": 0.01,
            "CO2": 0.00,
        }
        self.manual_toxin_sim = {
            "polysaccharide_lps": 0.20,
            "flagellin": 0.05,
        }

    def score_molecule(self, molecule_id: str, molecule_name: str,
                       smiles: str, molecular_weight: float = None) -> VaccinePotentialScore:
        """Enhanced scoring: boost vaccine candidates, penalise common metabolites."""
        rationale = []
        name_lower = molecule_name.lower()

        # ---- Boost for real vaccine candidates ----
        boost = 0.0
        if any(kw in name_lower for kw in ["lps", "lipopolysaccharide", "flagellin", "cpg", "polyic"]):
            boost += 0.6
            rationale.append("+0.6: Strong vaccine adjuvant pattern")
        if any(kw in name_lower for kw in ["virus", "capsid", "spike", "mannan", "beta_glucan", "vesicle"]):
            boost += 0.5
            rationale.append("+0.5: Structural vaccine component")
        if any(kw in name_lower for kw in ["toxoid"]):
            boost += 0.2
            rationale.append("+0.2: Toxoid")

        # ---- Penalty for common metabolites ----
        penalty = 0.0
        common = ["adenosine triphosphate", "atp", "glucose", "co2", "ethanol", "glycerol", "acetate", "lactate"]
        if name_lower in common:
            penalty = 0.7
            rationale.append("-0.7: Common metabolite (not vaccine)")

        # ---- Manual similarity scores ----
        sim_adj = self.manual_adjuvant_sim.get(molecule_id, 0.0)
        sim_tox = self.manual_toxin_sim.get(molecule_id, 0.0)
        rationale.append(f"Adjuvant similarity: {sim_adj:.3f}")
        rationale.append(f"Toxin similarity: {sim_tox:.3f}")

        # ---- Immunogenicity ----
        base_immuno = 0.3
        immunogenicity = max(0.0, min(1.0, base_immuno + boost - penalty))
        rationale.append(f"Immunogenicity: {immunogenicity:.3f}")

        # ---- Novelty (1 - similarity) ----
        novelty = 1.0 - sim_adj
        rationale.append(f"Novelty: {novelty:.3f}")

        # ---- Toxicity risk ----
        toxicity_risk = sim_tox

        # ---- Composite score (weights: immunogenicity 60, novelty 20, safety 10, adj 10) ----
        overall = (immunogenicity * 60 +
                   novelty * 20 +
                   (1 - toxicity_risk) * 10 +
                   sim_adj * 10)

        # ---- Recommendation ----
        if toxicity_risk > 0.7:
            recommendation = "excluded"
        elif overall >= 70:
            recommendation = "high"
        elif overall >= 35:
            recommendation = "medium"
        else:
            recommendation = "low"

        return VaccinePotentialScore(
            molecule_id=molecule_id,
            molecule_name=molecule_name,
            immunogenicity_score=immunogenicity,
            novelty_score=novelty,
            toxicity_risk=toxicity_risk,
            similarity_to_adjuvants=sim_adj,
            similarity_to_toxins=sim_tox,
            overall_vaccine_score=overall,
            recommendation=recommendation,
            rationale=rationale
        )

    def batch_score(self, molecules):
        """Score multiple molecules."""
        return [self.score_molecule(mid, name, smi, mw) for (mid, name, smi, mw) in molecules]

    # Dummy methods for compatibility (not used)
    def _similarity_to_known_set(self, *args, **kwargs): return 0.0
    def _tanimoto_similarity(self, *args, **kwargs): return 0.0
    def _string_similarity(self, *args, **kwargs): return 0.0
    def _estimate_immunogenicity(self, *args, **kwargs): return 0.3
    def _physicochemical_properties(self, *args, **kwargs): return 0.5, []


def create_default_scorer() -> VaccineScorer:
    """Factory function."""
    return VaccineScorer()