"""
MicroMol Combination Simulator - Core Engine
Real computational biology: metabolic networks, enzyme kinetics, combination screening.
"""

import json
import numpy as np
from typing import Dict, List, Tuple, Set, Optional
from dataclasses import dataclass, asdict
from collections import defaultdict
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class Metabolite:
    """Represents a small molecule or metabolite."""
    id: str
    name: str
    smiles: str  # Simplified Molecular Input Line Entry System
    formula: str
    mw: float  # Molecular weight
    can_be_secreted: bool = True
    source: str = "unknown"  # e.g., "nutrient", "product", "cofactor"


@dataclass
class Reaction:
    """Represents a biochemical reaction (enzyme-catalyzed)."""
    id: str
    name: str
    enzyme: str
    substrate_ids: List[str]
    product_ids: List[str]
    reversible: bool = False
    gpr: str = ""  # Gene-protein-reaction rule (e.g., "gene1 OR gene2")
    cofactors: Dict[str, float] = None  # e.g., {"ATP": 1, "NADH": 2}
    km_values: Dict[str, float] = None  # Michaelis constant for each substrate
    vmax: float = 1000.0  # Max velocity (micromol/min)
    threshold_for_activation: float = 0.01  # Min substrate conc for reaction to occur

    def __post_init__(self):
        if self.cofactors is None:
            self.cofactors = {}
        if self.km_values is None:
            self.km_values = {}


@dataclass
class Microbe:
    """Represents a microorganism with metabolic capabilities."""
    id: str
    name: str
    domain: str  # "bacteria", "virus", "fungus"
    reactions: Dict[str, Reaction]  # reaction_id -> Reaction
    native_metabolites: Set[str]  # Metabolites naturally produced by this microbe
    essential_nutrients: Set[str]  # What it MUST have to survive
    max_growth_rate: float = 0.5  # per hour
    
    def get_producible_metabolites(self) -> Set[str]:
        """All metabolites this microbe can produce (as products or cofactors)."""
        products = set()
        for rxn in self.reactions.values():
            products.update(rxn.product_ids)
        return products | self.native_metabolites

    def can_produce(self, metabolite_id: str) -> bool:
        """Can this microbe produce this metabolite?"""
        return metabolite_id in self.get_producible_metabolites()


class MetabolicNetwork:
    """Represents a complete metabolic network (single microbe or co-culture)."""
    
    def __init__(self, microbes: List[Microbe], shared_metabolites: Dict[str, Metabolite],
                 initial_conditions: Dict[str, float] = None):
        """
        Args:
            microbes: List of microbes in the culture
            shared_metabolites: Dict of all metabolites available
            initial_conditions: Dict mapping metabolite_id -> concentration (mM)
        """
        self.microbes = microbes
        self.metabolites = shared_metabolites
        self.concentrations = initial_conditions or {}
        self.time_steps = []
        self.concentration_history = {}
        
        # Initialize all metabolites to 0 if not specified
        for met_id in shared_metabolites:
            if met_id not in self.concentrations:
                self.concentrations[met_id] = 0.0

        # Biomass multipliers per microbe id (affects reaction fluxes)
        # Defaults to 1.0 for each supplied microbe
        self.biomass = {m.id: 1.0 for m in microbes}
    
    def simulate(self, duration: float = 24.0, dt: float = 0.1, 
                 nutrient_feed: Optional[Dict[str, float]] = None) -> Dict:
        """
        Simulate the metabolic network over time using simple Monod kinetics.
        
        Args:
            duration: Simulation time (hours)
            dt: Time step (hours)
            nutrient_feed: Dict of metabolite_id -> concentration_added_per_step (mM)
        
        Returns:
            Dict with final concentrations, history, and produced metabolites
        """
        if nutrient_feed is None:
            nutrient_feed = {}
        
        # Track which metabolites were present at start
        initial_metabolites = set(mid for mid, conc in self.concentrations.items() if conc > 0)
        
        num_steps = int(duration / dt)
        self.concentration_history = {mid: [] for mid in self.metabolites.keys()}
        self.time_steps = []
        
        for step in range(num_steps):
            current_time = step * dt
            self.time_steps.append(current_time)
            
            # Add nutrients
            for met_id, amount in nutrient_feed.items():
                self.concentrations[met_id] = self.concentrations.get(met_id, 0) + amount
            
            # Evaluate all active reactions in all microbes
            reaction_fluxes = self._calculate_reaction_fluxes()
            
            # Update metabolite concentrations based on reaction fluxes
            for rxn_id, flux in reaction_fluxes.items():
                # Find which reaction this is
                for microbe in self.microbes:
                    if rxn_id in microbe.reactions:
                        rxn = microbe.reactions[rxn_id]
                        # Decrease substrates
                        for substrate_id in rxn.substrate_ids:
                            self.concentrations[substrate_id] = max(
                                0, self.concentrations.get(substrate_id, 0) - flux * dt
                            )
                        # Increase products
                        for product_id in rxn.product_ids:
                            self.concentrations[product_id] = (
                                self.concentrations.get(product_id, 0) + flux * dt
                            )
                        break
            
            # Decaying: some metabolites are consumed by microbes
            for met_id in list(self.concentrations.keys()):
                if self.concentrations[met_id] > 0:
                    # Natural decay
                    self.concentrations[met_id] *= 0.99
            
            # Record history
            for met_id in self.metabolites.keys():
                self.concentration_history[met_id].append(
                    self.concentrations.get(met_id, 0)
                )
        
        # Identify newly produced metabolites
        final_metabolites = set(mid for mid, conc in self.concentrations.items() if conc > 1e-6)
        novel_metabolites = final_metabolites - initial_metabolites
        
        return {
            "final_concentrations": dict(self.concentrations),
            "novel_metabolites": list(novel_metabolites),
            "duration": duration,
            "time_steps": self.time_steps,
            "concentration_history": self.concentration_history
        }
    
    def _calculate_reaction_fluxes(self) -> Dict[str, float]:
        """
        Calculate reaction fluxes using Monod kinetics.
        
        Returns:
            Dict mapping reaction_id -> flux (micromol/min)
        """
        fluxes = {}
        
        for microbe in self.microbes:
            # biomass multiplier for this microbe
            bm = self.biomass.get(microbe.id, 1.0)
            for rxn_id, rxn in microbe.reactions.items():
                # Check if all substrates are present above threshold
                can_proceed = True
                for substrate_id in rxn.substrate_ids:
                    conc = self.concentrations.get(substrate_id, 0)
                    if conc < rxn.threshold_for_activation:
                        can_proceed = False
                        break
                
                if not can_proceed:
                    fluxes[rxn_id] = 0.0
                    continue
                
                # Monod kinetics: v = Vmax * [S] / (Km + [S])
                # For multiple substrates, use multiplicative form
                relative_velocity = 1.0
                for substrate_id in rxn.substrate_ids:
                    conc = self.concentrations.get(substrate_id, 0)
                    km = rxn.km_values.get(substrate_id, 0.1)
                    relative_velocity *= conc / (km + conc)
                
                flux = rxn.vmax * relative_velocity
                # scale flux by microbe biomass (initial inoculum or relative abundance)
                fluxes[rxn_id] = fluxes.get(rxn_id, 0.0) + flux * bm
        
        return fluxes


class CombinationSimulator:
    """Orchestrates simulation of microbe-molecule combinations."""
    
    def __init__(self, microbes: Dict[str, Microbe], metabolites: Dict[str, Metabolite]):
        """
        Args:
            microbes: Dict mapping microbe_id -> Microbe object
            metabolites: Dict mapping metabolite_id -> Metabolite object
        """
        self.microbes = microbes
        self.metabolites = metabolites
        self.results = []
    
    def simulate_monoculture(self, microbe_id: str, 
                            initial_nutrients: Dict[str, float],
                            duration: float = 24.0,
                            biomass_ratios: Optional[Dict[str, float]] = None) -> Dict:
        """
        Simulate a single microbe in isolation.
        
        Returns:
            Dict with metabolite concentrations and novel products
        """
        if microbe_id not in self.microbes:
            raise ValueError(f"Microbe {microbe_id} not found")
        
        microbe = self.microbes[microbe_id]
        network = MetabolicNetwork([microbe], self.metabolites, initial_nutrients.copy())
        # Apply biomass ratios if provided (override defaults)
        if biomass_ratios:
            for mid, val in biomass_ratios.items():
                if mid in network.biomass:
                    network.biomass[mid] = float(val)
        result = network.simulate(duration=duration, dt=0.1)
        result["microbes"] = [microbe_id]
        result["type"] = "monoculture"
        
        return result
    
    def simulate_coculture(self, microbe_ids: List[str],
                          initial_nutrients: Dict[str, float],
                          duration: float = 24.0,
                          allow_exchange: bool = True,
                          biomass_ratios: Optional[Dict[str, float]] = None) -> Dict:
        """
        Simulate a co-culture of multiple microbes.
        
        Args:
            microbe_ids: List of microbe IDs
            initial_nutrients: Initial metabolite concentrations
            duration: Simulation time (hours)
            allow_exchange: Whether microbes can exchange metabolites
        
        Returns:
            Dict with simulation results
        """
        microbes = [self.microbes[mid] for mid in microbe_ids if mid in self.microbes]
        
        if len(microbes) == 0:
            raise ValueError("No valid microbes specified")
        
        network = MetabolicNetwork(microbes, self.metabolites, initial_nutrients.copy())
        # If biomass ratios passed, apply them to the network. Missing entries keep default 1.0
        if biomass_ratios:
            for mid, val in biomass_ratios.items():
                if mid in network.biomass:
                    network.biomass[mid] = float(val)
        result = network.simulate(duration=duration, dt=0.1)
        result["microbes"] = microbe_ids
        result["type"] = "coculture"
        result["allow_exchange"] = allow_exchange
        
        return result
    
    def simulate_with_substrate(self, microbe_ids: List[str], 
                               substrate_id: str,
                               substrate_conc: float = 10.0,
                               duration: float = 24.0) -> Dict:
        """
        Simulate with a specific substrate and track products.
        
        Returns:
            Result dict with produced metabolites
        """
        initial = {substrate_id: substrate_conc}
        return self.simulate_coculture(microbe_ids, initial, duration)


# ==============================================================================
# PRESET DATABASES: Simplified but realistic metabolic models
# ==============================================================================

def create_ecoli_reactions() -> Dict[str, Reaction]:
    """
    E. coli core metabolic reactions (simplified).
    Real model would use BiGG or ModelSEED.
    """
    return {
        "r1_glycolysis_fwd": Reaction(
            id="r1_glycolysis_fwd",
            name="Glucose → Pyruvate (Glycolysis)",
            enzyme="Hexokinase/PFK complex",
            substrate_ids=["glucose", "ATP"],
            product_ids=["pyruvate", "ADP", "NADH"],
            cofactors={"ATP": -2, "ADP": 2, "NADH": 2},
            km_values={"glucose": 0.5, "ATP": 0.1},
            vmax=100.0
        ),
        "r2_acetate_formation": Reaction(
            id="r2_acetate_formation",
            name="Pyruvate → Acetate (overflow metabolism)",
            enzyme="Pyruvate dehydrogenase complex",
            substrate_ids=["pyruvate"],
            product_ids=["acetate", "NADH"],
            km_values={"pyruvate": 0.2},
            vmax=50.0
        ),
        "r3_lactate_from_pyruvate": Reaction(
            id="r3_lactate_from_pyruvate",
            name="Pyruvate → Lactate (anaerobic)",
            enzyme="Lactate dehydrogenase",
            substrate_ids=["pyruvate", "NADH"],
            product_ids=["lactate"],
            km_values={"pyruvate": 0.3, "NADH": 0.05},
            vmax=40.0
        ),
        "r4_alanine_synthesis": Reaction(
            id="r4_alanine_synthesis",
            name="Pyruvate → Alanine (transamination)",
            enzyme="Alanine aminotransferase",
            substrate_ids=["pyruvate", "glutamate"],
            product_ids=["alanine", "alpha_ketoglutarate"],
            km_values={"pyruvate": 0.2, "glutamate": 0.1},
            vmax=30.0
        ),
        "r5_polysaccharide": Reaction(
            id="r5_polysaccharide",
            name="Glucose → Polysaccharide (biofilm)",
            enzyme="Glycosyltransferase",
            substrate_ids=["glucose"],
            product_ids=["polysaccharide_lps"],
            km_values={"glucose": 0.4},
            vmax=20.0
        ),
    }


def create_yeast_reactions() -> Dict[str, Reaction]:
    """
    Saccharomyces cerevisiae fermentation reactions (simplified).
    """
    return {
        "y1_glucose_fermentation": Reaction(
            id="y1_glucose_fermentation",
            name="Glucose → Ethanol (fermentation)",
            enzyme="Pyruvate decarboxylase + Alcohol dehydrogenase",
            substrate_ids=["glucose"],
            product_ids=["ethanol", "CO2"],
            km_values={"glucose": 0.3},
            vmax=80.0
        ),
        "y2_glucose_respiration": Reaction(
            id="y2_glucose_respiration",
            name="Glucose → CO2 (aerobic respiration)",
            enzyme="Citric acid cycle",
            substrate_ids=["glucose"],
            product_ids=["CO2", "ATP"],
            km_values={"glucose": 0.2},
            vmax=60.0
        ),
        "y3_glycerol_production": Reaction(
            id="y3_glycerol_production",
            name="Glucose → Glycerol (osmolyte)",
            enzyme="Glycerol-3-phosphate dehydrogenase",
            substrate_ids=["glucose"],
            product_ids=["glycerol"],
            km_values={"glucose": 0.4},
            vmax=25.0
        ),
        "y4_polysaccharide": Reaction(
            id="y4_polysaccharide",
            name="Glucose → Mannan (cell wall)",
            enzyme="Mannosyltransferase",
            substrate_ids=["glucose"],
            product_ids=["mannan_polysaccharide"],
            km_values={"glucose": 0.35},
            vmax=15.0
        ),
    }


def create_virus_reactions() -> Dict[str, Reaction]:
    """
    Simplified viral reaction set (depends on host nucleotides).
    Viruses don't have own metabolism but hijack host machinery.
    """
    return {
        "v1_capsid_assembly": Reaction(
            id="v1_capsid_assembly",
            name="Host proteins → Viral capsid",
            enzyme="Viral protease + assembly factors",
            substrate_ids=["host_proteins"],
            product_ids=["viral_capsid"],
            km_values={"host_proteins": 0.2},
            vmax=10.0
        ),
        "v2_rna_replication": Reaction(
            id="v2_rna_replication",
            name="Host nucleotides → Viral RNA",
            enzyme="Viral RNA polymerase",
            substrate_ids=["host_nucleotides"],
            product_ids=["viral_rna"],
            km_values={"host_nucleotides": 0.15},
            vmax=15.0
        ),
        "v3_surface_protein_display": Reaction(
            id="v3_surface_protein_display",
            name="Viral proteins on capsid surface",
            enzyme="Viral assembly",
            substrate_ids=["viral_spike_protein"],
            product_ids=["viral_particle_with_spike"],
            km_values={"viral_spike_protein": 0.1},
            vmax=8.0
        ),
    }


def create_metabolite_registry() -> Dict[str, Metabolite]:
    """Complete metabolite database with SMILES and properties."""
    return {
        # Nutrients
        "glucose": Metabolite(
            id="glucose",
            name="Glucose",
            smiles="C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O",
            formula="C6H12O6",
            mw=180.16,
            source="nutrient"
        ),
        "acetate": Metabolite(
            id="acetate",
            name="Acetic acid",
            smiles="CC(=O)O",
            formula="C2H4O2",
            mw=60.05,
            source="nutrient"
        ),
        # Core metabolites
        "pyruvate": Metabolite(
            id="pyruvate",
            name="Pyruvate",
            smiles="CC(=O)[O-]",
            formula="C3H3O3",
            mw=87.03,
            source="metabolite"
        ),
        "lactate": Metabolite(
            id="lactate",
            name="Lactic acid",
            smiles="CC(O)C(=O)O",
            formula="C3H6O3",
            mw=90.08,
            source="metabolite"
        ),
        "alanine": Metabolite(
            id="alanine",
            name="Alanine",
            smiles="CC(N)C(=O)O",
            formula="C3H7NO2",
            mw=89.09,
            source="metabolite"
        ),
        "glutamate": Metabolite(
            id="glutamate",
            name="Glutamic acid",
            smiles="C(CC(=O)O)[C@@H](C(=O)O)N",
            formula="C5H9NO4",
            mw=147.13,
            source="metabolite"
        ),
        "alpha_ketoglutarate": Metabolite(
            id="alpha_ketoglutarate",
            name="α-Ketoglutarate",
            smiles="C(CC(=O)[O-])C(=O)C(=O)[O-]",
            formula="C5H6O5",
            mw=146.10,
            source="metabolite"
        ),
        # Fermentation products
        "ethanol": Metabolite(
            id="ethanol",
            name="Ethanol",
            smiles="CCO",
            formula="C2H6O",
            mw=46.07,
            source="metabolite"
        ),
        "glycerol": Metabolite(
            id="glycerol",
            name="Glycerol",
            smiles="C(C(CO)O)O",
            formula="C3H8O3",
            mw=92.09,
            source="metabolite"
        ),
        # Cofactors
        "ATP": Metabolite(
            id="ATP",
            name="Adenosine triphosphate",
            smiles="C1=NC(=C2C(=N1)N(C=N2)[C@@H]3[C@H]([C@H]([C@@H](O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N",
            formula="C10H16N5O13P3",
            mw=507.18,
            source="cofactor"
        ),
        "NADH": Metabolite(
            id="NADH",
            name="NADH",
            smiles="CC(=O)N[C@@H]1[C@H]([C@H]([C@@H](O1)COP(=O)([O-])OP(=O)([O-])OC[C@@H]2[C@H]([C@@H]([C@H](O2)N3C=NC4=C3N=CN=C4N)O)OP(=O)([O-])[O-])[O-])O",
            formula="C21H29N7O14P2",
            mw=665.42,
            source="cofactor"
        ),
        "ADP": Metabolite(
            id="ADP",
            name="Adenosine diphosphate",
            smiles="C1=NC(=C2C(=N1)N(C=N2)[C@@H]3[C@H]([C@H]([C@@H](O3)COP(=O)(O)OP(=O)(O)O)O)O)N",
            formula="C10H15N5O10P2",
            mw=427.20,
            source="cofactor"
        ),
        "CO2": Metabolite(
            id="CO2",
            name="Carbon dioxide",
            smiles="C(=O)=O",
            formula="CO2",
            mw=44.01,
            source="gas"
        ),
        # Bacterial metabolites
        "polysaccharide_lps": Metabolite(
            id="polysaccharide_lps",
            name="LPS (Lipopolysaccharide)",
            smiles="CC(C)CC(C)(C)C(=O)O[C@@H]1[C@@H](O)[C@H](O)[C@@H](C)O[C@H]1OC[C@@]2(O)[C@H](O)[C@@H](O)[C@H](CO)O2",
            formula="C20H36O13",  # Simplified; real LPS is much larger
            mw=468.21,
            source="metabolite"
        ),
        # Yeast metabolites
        "mannan_polysaccharide": Metabolite(
            id="mannan_polysaccharide",
            name="Mannan",
            smiles="C([C@@H]1[C@H]([C@H]([C@H](C(O1)O[C@@H]2[C@H]([C@@H]([C@H](C(O2)O)O)O)O)O)O)O)O",
            formula="C12H22O11",  # Simplified
            mw=342.30,
            source="metabolite"
        ),
        # Viral-related
        "host_proteins": Metabolite(
            id="host_proteins",
            name="Host proteins (substrate)",
            smiles="C(N)CC(=O)O",  # Simplified as amino acid
            formula="C3H7NO2",
            mw=89.09,
            source="metabolite"
        ),
        "host_nucleotides": Metabolite(
            id="host_nucleotides",
            name="Host nucleotides",
            smiles="C1=C[N+](=CC=C1)[C@@H]2[C@H]([C@H]([C@@H](O2)COP(=O)(O)O)O)O",  # Simplified
            formula="C10H14N4O7P",
            mw=305.21,
            source="metabolite"
        ),
        "viral_rna": Metabolite(
            id="viral_rna",
            name="Viral RNA",
            smiles="C(C(C(C(CO)O)O)O)O",  # Simplified as ribose polymer
            formula="C5H10O5",
            mw=150.13,
            source="metabolite"
        ),
        "viral_capsid": Metabolite(
            id="viral_capsid",
            name="Viral capsid protein assembly",
            smiles="C(CC(N)C(=O)O)(C(=O)O)N",  # Simplified
            formula="C5H10N2O4",
            mw=162.15,
            source="metabolite"
        ),
        "viral_spike_protein": Metabolite(
            id="viral_spike_protein",
            name="Viral spike protein",
            smiles="C(C(C(C(CO)O)O)O)O",  # Simplified
            formula="C5H12O6",
            mw=152.15,
            source="metabolite"
        ),
        "viral_particle_with_spike": Metabolite(
            id="viral_particle_with_spike",
            name="Viral particle with spike protein",
            smiles="C(C(C(C(CO)O)O)O)O",  # Very simplified
            formula="C10H20O12",
            mw=300.30,
            source="metabolite"
        ),
    }


def create_default_microbes() -> Dict[str, Microbe]:
    """Create default microbe registry."""
    return {
        "ecoli": Microbe(
            id="ecoli",
            name="Escherichia coli",
            domain="bacteria",
            reactions=create_ecoli_reactions(),
            native_metabolites={"acetate", "lactate", "alanine"},
            essential_nutrients={"glucose"},
            max_growth_rate=0.7
        ),
        "yeast": Microbe(
            id="yeast",
            name="Saccharomyces cerevisiae",
            domain="fungus",
            reactions=create_yeast_reactions(),
            native_metabolites={"ethanol", "glycerol"},
            essential_nutrients={"glucose"},
            max_growth_rate=0.4
        ),
        "flu_virus": Microbe(
            id="flu_virus",
            name="Influenza virus",
            domain="virus",
            reactions=create_virus_reactions(),
            native_metabolites={"viral_rna", "viral_capsid"},
            essential_nutrients={"host_proteins", "host_nucleotides"},
            max_growth_rate=0.3
        ),
        # Additional custom microbes that produce vaccine-relevant molecules
        "bacillus": Microbe(
            id="bacillus",
            name="Bacillus subtilis",
            domain="bacteria",
            reactions={
                **create_ecoli_reactions(),
                "flagellin_synthesis": Reaction(
                    id="flagellin_synthesis",
                    name="Flagellin synthesis",
                    enzyme="Flagellin synthase",
                    substrate_ids=["alanine", "ATP"],
                    product_ids=["flagellin"],
                    km_values={"alanine": 0.5, "ATP": 0.2},
                    vmax=20.0
                ),
                "teichoic_acid_synthesis": Reaction(
                    id="teichoic_acid_synthesis",
                    name="Teichoic acid synthesis",
                    enzyme="TagF",
                    substrate_ids=["glycerol", "ATP"],
                    product_ids=["lipoteichoic_acid"],
                    km_values={"glycerol": 0.3, "ATP": 0.1},
                    vmax=15.0
                )
            },
            native_metabolites={"flagellin", "lipoteichoic_acid"},
            essential_nutrients={"glucose"},
            max_growth_rate=0.6
        ),

        "salmonella": Microbe(
            id="salmonella",
            name="Salmonella typhimurium",
            domain="bacteria",
            reactions={
                **create_ecoli_reactions(),
                "lps_synthesis": Reaction(
                    id="lps_synthesis",
                    name="LPS synthesis",
                    enzyme="Lpx acyltransferase",
                    substrate_ids=["glucose", "ATP", "acetate"],
                    product_ids=["polysaccharide_lps"],
                    km_values={"glucose": 0.2, "ATP": 0.1, "acetate": 0.05},
                    vmax=25.0
                ),
                "flagellin_synthesis": Reaction(
                    id="flagellin_synthesis",
                    name="Flagellin synthesis",
                    enzyme="Flagellin synthase",
                    substrate_ids=["alanine", "ATP"],
                    product_ids=["flagellin"],
                    km_values={"alanine": 0.5, "ATP": 0.2},
                    vmax=18.0
                )
            },
            native_metabolites={"polysaccharide_lps", "flagellin"},
            essential_nutrients={"glucose"},
            max_growth_rate=0.7
        ),

        "lactobacillus": Microbe(
            id="lactobacillus",
            name="Lactobacillus plantarum",
            domain="bacteria",
            reactions={
                **create_yeast_reactions(),
                "mannan_synthesis": Reaction(
                    id="mannan_synthesis",
                    name="Mannan synthesis",
                    enzyme="Mannosyltransferase",
                    substrate_ids=["glucose"],
                    product_ids=["mannan_polysaccharide"],
                    km_values={"glucose": 0.2},
                    vmax=30.0
                ),
                "peptidoglycan_synthesis": Reaction(
                    id="peptidoglycan_synthesis",
                    name="Peptidoglycan synthesis",
                    enzyme="Mur ligase",
                    substrate_ids=["glucose", "alanine"],
                    product_ids=["peptidoglycan"],
                    km_values={"glucose": 0.3, "alanine": 0.4},
                    vmax=20.0
                )
            },
            native_metabolites={"mannan_polysaccharide", "peptidoglycan"},
            essential_nutrients={"glucose"},
            max_growth_rate=0.5
        ),
        # Ten additional simple custom microbes to expand search space
        "bacteroides": Microbe(
            id="bacteroides",
            name="Bacteroides thetaiotaomicron",
            domain="bacteria",
            reactions={
                **create_ecoli_reactions(),
                "omv_release": Reaction(
                    id="omv_release",
                    name="Outer membrane vesicle release",
                    enzyme="OMV machinery",
                    substrate_ids=["acetate"],
                    product_ids=["outer_membrane_vesicle"],
                    km_values={"acetate": 0.1},
                    vmax=10.0
                )
            },
            native_metabolites={"outer_membrane_vesicle"},
            essential_nutrients={"glucose"},
            max_growth_rate=0.45
        ),

        "bifidobacterium": Microbe(
            id="bifidobacterium",
            name="Bifidobacterium longum",
            domain="bacteria",
            reactions={
                **create_yeast_reactions(),
                "exopolysaccharide": Reaction(
                    id="exopolysaccharide",
                    name="Exopolysaccharide production",
                    enzyme="EPS synthase",
                    substrate_ids=["glucose"],
                    product_ids=["exopolysaccharide"],
                    km_values={"glucose": 0.25},
                    vmax=18.0
                )
            },
            native_metabolites={"exopolysaccharide"},
            essential_nutrients={"glucose"},
            max_growth_rate=0.35
        ),

        "clostridium": Microbe(
            id="clostridium",
            name="Clostridium difficile",
            domain="bacteria",
            reactions={
                **create_ecoli_reactions(),
                "butyrate_production": Reaction(
                    id="butyrate_production",
                    name="Butyrate production",
                    enzyme="Butyrate kinase",
                    substrate_ids=["pyruvate"],
                    product_ids=["butyrate"],
                    km_values={"pyruvate": 0.2},
                    vmax=22.0
                )
            },
            native_metabolites={"butyrate"},
            essential_nutrients={"glucose"},
            max_growth_rate=0.4
        ),

        "streptococcus": Microbe(
            id="streptococcus",
            name="Streptococcus pneumoniae",
            domain="bacteria",
            reactions={
                **create_ecoli_reactions(),
                "capsular_polysaccharide": Reaction(
                    id="capsular_polysaccharide",
                    name="Capsular polysaccharide synthesis",
                    enzyme="Cps synthase",
                    substrate_ids=["glucose"],
                    product_ids=["capsular_polysaccharide"],
                    km_values={"glucose": 0.3},
                    vmax=16.0
                )
            },
            native_metabolites={"capsular_polysaccharide"},
            essential_nutrients={"glucose"},
            max_growth_rate=0.5
        ),

        "enterococcus": Microbe(
            id="enterococcus",
            name="Enterococcus faecalis",
            domain="bacteria",
            reactions={
                **create_ecoli_reactions(),
                "lipoprotein_release": Reaction(
                    id="lipoprotein_release",
                    name="Lipoprotein release",
                    enzyme="Lipoprotein synthase",
                    substrate_ids=["alanine"],
                    product_ids=["bacterial_lipoprotein"],
                    km_values={"alanine": 0.2},
                    vmax=12.0
                )
            },
            native_metabolites={"bacterial_lipoprotein"},
            essential_nutrients={"glucose"},
            max_growth_rate=0.45
        ),

        "staphylococcus": Microbe(
            id="staphylococcus",
            name="Staphylococcus aureus",
            domain="bacteria",
            reactions={
                **create_ecoli_reactions(),
                "teichoic_release": Reaction(
                    id="teichoic_release",
                    name="Teichoic acid shedding",
                    enzyme="Tag synthase",
                    substrate_ids=["glycerol"],
                    product_ids=["teichoic_acid"],
                    km_values={"glycerol": 0.3},
                    vmax=14.0
                )
            },
            native_metabolites={"teichoic_acid"},
            essential_nutrients={"glucose"},
            max_growth_rate=0.55
        ),

        "pseudomonas": Microbe(
            id="pseudomonas",
            name="Pseudomonas aeruginosa",
            domain="bacteria",
            reactions={
                **create_ecoli_reactions(),
                "alginate_production": Reaction(
                    id="alginate_production",
                    name="Alginate production",
                    enzyme="Alginate synthase",
                    substrate_ids=["glucose"],
                    product_ids=["alginate"],
                    km_values={"glucose": 0.35},
                    vmax=18.0
                )
            },
            native_metabolites={"alginate"},
            essential_nutrients={"glucose"},
            max_growth_rate=0.6
        ),

        "klebsiella": Microbe(
            id="klebsiella",
            name="Klebsiella pneumoniae",
            domain="bacteria",
            reactions={
                **create_ecoli_reactions(),
                "capsule_synthesis": Reaction(
                    id="capsule_synthesis",
                    name="Capsule polysaccharide synthesis",
                    enzyme="Wzy polymerase",
                    substrate_ids=["glucose"],
                    product_ids=["capsule_polysaccharide"],
                    km_values={"glucose": 0.3},
                    vmax=20.0
                )
            },
            native_metabolites={"capsule_polysaccharide"},
            essential_nutrients={"glucose"},
            max_growth_rate=0.5
        ),

        "faecalibacterium": Microbe(
            id="faecalibacterium",
            name="Faecalibacterium prausnitzii",
            domain="bacteria",
            reactions={
                **create_ecoli_reactions(),
                "butyrate_from_acetate": Reaction(
                    id="butyrate_from_acetate",
                    name="Butyrate from acetate",
                    enzyme="Butyryl-CoA:acetate CoA-transferase",
                    substrate_ids=["acetate"],
                    product_ids=["butyrate"],
                    km_values={"acetate": 0.1},
                    vmax=25.0
                )
            },
            native_metabolites={"butyrate"},
            essential_nutrients={"glucose"},
            max_growth_rate=0.4
        ),
        # Five more microbes to expand search
        "mycobacterium": Microbe(
            id="mycobacterium",
            name="Mycobacterium smegmatis",
            domain="bacteria",
            reactions={
                **create_ecoli_reactions(),
                "mycolic_acid_synthesis": Reaction(
                    id="mycolic_acid_synthesis",
                    name="Mycolic acid synthesis",
                    enzyme="Mycolyl transferase",
                    substrate_ids=["fatty_acids"],
                    product_ids=["mycolic_acid"],
                    km_values={"fatty_acids": 0.2},
                    vmax=12.0
                )
            },
            native_metabolites={"mycolic_acid"},
            essential_nutrients={"glucose"},
            max_growth_rate=0.3
        ),

        "bordetella": Microbe(
            id="bordetella",
            name="Bordetella pertussis",
            domain="bacteria",
            reactions={
                **create_ecoli_reactions(),
                "pertussis_toxin_component": Reaction(
                    id="pertussis_toxin_component",
                    name="Pertussis toxin subunit synthesis",
                    enzyme="Toxin synth",
                    substrate_ids=["alanine"],
                    product_ids=["pertussis_subunit"],
                    km_values={"alanine": 0.2},
                    vmax=8.0
                )
            },
            native_metabolites={"pertussis_subunit"},
            essential_nutrients={"glucose"},
            max_growth_rate=0.35
        ),

        "vibrio": Microbe(
            id="vibrio",
            name="Vibrio cholerae",
            domain="bacteria",
            reactions={
                **create_ecoli_reactions(),
                "cholera_toxin_component": Reaction(
                    id="cholera_toxin_component",
                    name="Cholera toxin subunit synthesis",
                    enzyme="Ctx synth",
                    substrate_ids=["alanine"],
                    product_ids=["cholera_subunit"],
                    km_values={"alanine": 0.2},
                    vmax=10.0
                )
            },
            native_metabolites={"cholera_subunit"},
            essential_nutrients={"glucose"},
            max_growth_rate=0.5
        ),

        "helicobacter": Microbe(
            id="helicobacter",
            name="Helicobacter pylori",
            domain="bacteria",
            reactions={
                **create_ecoli_reactions(),
                "urease_activity": Reaction(
                    id="urease_activity",
                    name="Urease activity",
                    enzyme="Urease",
                    substrate_ids=["urea"],
                    product_ids=["ammonia"],
                    km_values={"urea": 0.1},
                    vmax=20.0
                )
            },
            native_metabolites={"ammonia"},
            essential_nutrients={"glucose"},
            max_growth_rate=0.4
        ),

        "legionella": Microbe(
            id="legionella",
            name="Legionella pneumophila",
            domain="bacteria",
            reactions={
                **create_ecoli_reactions(),
                "dot_secretion": Reaction(
                    id="dot_secretion",
                    name="Secretion system effector",
                    enzyme="Dot/Icm",
                    substrate_ids=["host_proteins"],
                    product_ids=["effector_protein"],
                    km_values={"host_proteins": 0.2},
                    vmax=6.0
                )
            },
            native_metabolites={"effector_protein"},
            essential_nutrients={"glucose"},
            max_growth_rate=0.35
        ),
    }
