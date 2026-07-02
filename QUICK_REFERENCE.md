# MicroMol Simulator - Quick Reference & Examples

## COMMAND CHEAT SHEET

### Basic Commands

```bash
# Run full simulation (all microbes, 24h)
python run_simulator.py

# Two specific microbes only
python run_simulator.py --microbes ecoli yeast

# Longer simulation (48 hours)
python run_simulator.py --duration 48

# Different substrate (acetate at 5 mM)
python run_simulator.py --substrate acetate --concentration 5

# Use only 2 parallel workers
python run_simulator.py --workers 2

# No parallelization (serial mode, useful for debugging)
python run_simulator.py --serial

# Save to custom directory
python run_simulator.py --output-dir ~/experiments/exp001

# Skip monocultures, only combinations
python run_simulator.py --skip-mono

# All options together
python run_simulator.py \
  --microbes ecoli yeast bacillus \
  --duration 36 \
  --substrate acetate \
  --concentration 8.0 \
  --workers 4 \
  --output-dir ~/experiments/vaccine_screen_001
```

---

## FILE STRUCTURE REFERENCE

```
your_project/
├── micromol_core.py           # Engine: Monod kinetics, metabolic network
├── micromol_scorer.py         # Scoring: vaccine potential, chemical similarity
├── run_simulator.py           # Main executable: orchestration, I/O
├── SETUP_GUIDE.md            # Installation & usage (this folder)
├── TECHNICAL_ROADMAP.md      # Advanced features & next steps
├── micromol_results/         # Output directory (created automatically)
│   ├── micromol_results.csv          # Summary table
│   ├── micromol_detailed.json        # Full scoring details
│   └── micromol_summary.txt          # Human-readable report
└── models/                   # (Optional) BiGG SBML files
    ├── e_coli_core.xml
    ├── yeast_model.xml
    └── ...
```

---

## UNDERSTANDING THE OUTPUT FILES

### `micromol_results.csv`

```csv
Combination,Novel Metabolites Count,Novel Metabolites,Best Vaccine Score,Best Recommendation
ecoli + yeast,6,"glycerol, polysaccharide_lps, CO2, ethanol, ATP, mannan_polysaccharide",26.2,low
ecoli + flu_virus,3,"viral_rna, viral_capsid, viral_particle_with_spike",35.1,medium
yeast + flu_virus,2,"ethanol, CO2",12.5,low
```

**How to interpret:**
- **Novel Metabolites Count**: More ≠ better. Look for fewer, higher-quality candidates.
- **Best Vaccine Score**: 0-100 scale. >60 = "high" recommendation, 35-60 = "medium", <35 = "low"
- **Best Recommendation**: See scoring breakdown in JSON for details.

### `micromol_detailed.json`

```json
[
  {
    "combination": "ecoli + yeast",
    "microbe_ids": ["ecoli", "yeast"],
    "novel_metabolites": ["glycerol", "ATP", ...],
    "novel_metabolite_count": 6,
    "scores": [
      {
        "molecule_id": "glycerol",
        "molecule_name": "Glycerol",
        "immunogenicity_score": 0.25,
        "novelty_score": 0.45,
        "toxicity_risk": 0.08,
        "similarity_to_adjuvants": 0.12,
        "similarity_to_toxins": 0.02,
        "overall_vaccine_score": 26.23,
        "recommendation": "low",
        "rationale": [
          "Similarity to known adjuvants: 0.120",
          "Similarity to known toxins: 0.020",
          "Estimated immunogenicity: 0.250",
          "Novelty score (1 - avg similarity): 0.450",
          ...
        ]
      },
      ...
    ],
    "final_concentrations": {
      "glucose": 2.45,
      "pyruvate": 0.89,
      "ethanol": 5.32,
      ...
    }
  }
]
```

**Use `micromol_detailed.json` to:**
- Understand why a molecule got a certain score
- Check if scoring logic makes sense for your domain
- Aggregate results programmatically

### `micromol_summary.txt`

Human-readable report. Example:

```
======================================================================
MICROMOL COMBINATION SIMULATOR - SUMMARY REPORT
======================================================================

Timestamp: 2026-05-23 16:54:47

MONOCULTURE RESULTS (3 runs):
----------------------------------------------------------------------
  ecoli
    Novel metabolites: 2
  yeast
    Novel metabolites: 1
  flu_virus
    Novel metabolites: 0

COMBINATION RESULTS (3 combinations):
----------------------------------------------------------------------
TOP 5 CANDIDATES (by vaccine potential):

  1. ecoli + yeast
     Novel metabolites: 6
     Best candidate: Glycerol
     Vaccine score: 26.2/100
     Recommendation: low

  2. ecoli + flu_virus
     Novel metabolites: 3
     Best candidate: Viral RNA
     Vaccine score: 35.1/100
     Recommendation: medium

  3. yeast + flu_virus
     Novel metabolites: 2
     Best candidate: Ethanol
     Vaccine score: 12.5/100
     Recommendation: low

======================================================================
END OF REPORT
```

---

## EXAMPLE USE CASES

### Use Case 1: Screen for Novel COVID-19 Adjuvants

**Goal:** Find new molecules that SARS-CoV-2 produces in co-culture with human immune cells.

**Approach:**

```bash
# 1. Edit micromol_core.py: add human immune cell models
#    (dendritic cells, macrophages, B cells)

# 2. Create new adjuvant database in micromol_scorer.py
#    (known COVID vaccine adjuvants: AS01b, MF59, etc.)

# 3. Run screening
python run_simulator.py \
  --microbes sars_cov2 human_dendritic_cell macrophage \
  --duration 72 \
  --workers 8 \
  --output-dir covid_adjuvant_screen

# 4. Inspect results
cat covid_adjuvant_screen/micromol_summary.txt
# Look for: "high" or "medium" recommendations

# 5. Validate top 3 candidates in lab
#    (TLR activation assays, cytokine production, immunization)
```

**Expected outcome:**
- 3-5 novel adjuvant candidates
- Score correlation with wet-lab immunogenicity > 0.7

---

### Use Case 2: Optimize Probiotic Combination for Gut Health

**Goal:** Find which combination of gut bacteria produces the most immunoprotective metabolites.

**Approach:**

```bash
# 1. Load real 16S data from healthy humans
#    (download AGORA microbiome models: https://vmh.uni.lu/)

# 2. Screen common species combinations
python run_simulator.py \
  --microbes bacteroides_fragilis faecalibacterium_prausnitzii \
                roseburia_faecalis akkermansia_muciniphila \
  --duration 48 \
  --concentration 20 \
  --output-dir probiotic_screen

# 3. Check for short-chain fatty acids (SCFAs)
#    (butyrate, propionate → immune tolerance)
#    These should appear in "Novel Metabolites"

# 4. Rank combinations by SCFA production
#    (Manually: sum SCFA concentrations in final_concentrations)

# 5. Formulate probiotic supplement
#    (Capsule with top-ranked combination)
```

**Expected outcome:**
- Probiotic formulation producing 2-3 SCFAs per combination
- Better than random mixtures

---

### Use Case 3: Engineer Microbe for Heterologous Protein Production

**Goal:** Design a strain that produces a therapeutic protein in co-culture.

**Approach:**

```bash
# 1. Create engineered variant in micromol_core.py
microbes["ecoli_gfp"] = Microbe(
    id="ecoli_gfp",
    name="E. coli with GFP gene",
    reactions={
        **create_ecoli_reactions(),
        "gfp_synthesis": Reaction(
            id="gfp_synthesis",
            name="GFP protein synthesis",
            substrate_ids=["alanine", "ATP"],  # Simplified
            product_ids=["gfp_protein"],
            vmax=50
        ),
    },
    native_metabolites={"gfp_protein"}
)

# 2. Run with nutrient supplier
python run_simulator.py \
  --microbes ecoli_gfp yeast \
  --duration 48 \
  --substrate glucose \
  --concentration 30 \
  --output-dir protein_expression

# 3. Check if "gfp_protein" appears in novel_metabolites
# 4. Score it (should be "high" if engineered correctly)
```

---

### Use Case 4: Predict Antibiotic Resistance from Metabolite Exchange

**Goal:** Identify metabolites that confer antibiotic resistance in polymicrobial biofilms.

**Approach:**

```python
# Edit micromol_core.py: add resistance mechanism

"resistant_ecoli": Microbe(
    id="resistant_ecoli",
    name="Antibiotic-resistant E. coli",
    reactions={
        **create_ecoli_reactions(),
        "resistance_enzyme_production": Reaction(
            id="resistance_enzyme_production",
            name="Beta-lactamase synthesis",
            substrate_ids=["alanine"],
            product_ids=["beta_lactamase"],
            vmax=30
        ),
    }
)

# Run
python run_simulator.py \
  --microbes resistant_ecoli sensitive_pseudomonas \
  --duration 36 \
  --output-dir resistance_transfer
```

**Hypothesis:** In co-culture, resistant strain "feeds" the sensitive one with β-lactamase, conferring resistance via horizontal gene transfer.

---

### Use Case 5: High-Throughput Discovery (1,000+ combinations)

**Goal:** Screen all possible 3-way combinations of 20 microbes (~1,000 simulations).

**Approach:**

```bash
# Edit run_simulator.py:
# Change config.combination_sizes = [3]  # Triplets instead of pairs

# Prepare: download 20 BiGG models and place in models/

# Run on cluster
python run_simulator.py \
  --microbes ecoli yeast bacillus salmonella staph_aureus \
             listeria streptococcus pseudomonas acinetobacter \
             clostridium bacteroides prevotella faecalibacterium \
             roseburia akkermansia bifidobacterium lactobacillus \
             enterococcus \
  --duration 48 \
  --workers 32 \
  --output-dir discovery_1000_combos

# Results: ~10-50 "high" recommendation candidates
# Expected: some truly novel molecules never synthesized before
```

---

## CODE SNIPPETS FOR COMMON TASKS

### Add a New Microbe

```python
# In micromol_core.py, add to create_default_microbes():

"salmonella": Microbe(
    id="salmonella",
    name="Salmonella typhimurium",
    domain="bacteria",
    reactions=create_salmonella_reactions(),
    native_metabolites={"acetate", "lactate"},
    essential_nutrients={"glucose"},
    max_growth_rate=0.6
)

def create_salmonella_reactions() -> Dict[str, Reaction]:
    return {
        "s1_glucose_fermentation": Reaction(
            id="s1_glucose_fermentation",
            name="Glucose → Lactate (Salmonella)",
            enzyme="Lactate dehydrogenase",
            substrate_ids=["glucose"],
            product_ids=["lactate"],
            km_values={"glucose": 0.2},
            vmax=60.0
        ),
        # Add more Salmonella-specific reactions
    }
```

Then run:
```bash
python run_simulator.py --microbes ecoli salmonella
```

### Query Results Programmatically

```python
import json
import pandas as pd

# Load detailed JSON results
with open("micromol_results/micromol_detailed.json") as f:
    detailed = json.load(f)

# Find all high-scoring molecules
high_scores = []
for combo in detailed:
    for score in combo.get("scores", []):
        if score["recommendation"] == "high":
            high_scores.append({
                "combination": combo["combination"],
                "molecule": score["molecule_name"],
                "score": score["overall_vaccine_score"],
                "immunogenicity": score["immunogenicity_score"],
                "novelty": score["novelty_score"],
            })

# Convert to DataFrame
df = pd.DataFrame(high_scores)
df = df.sort_values("score", ascending=False)

# Export top 10
df.head(10).to_csv("top_candidates.csv", index=False)
print(df.head(10))
```

### Custom Scoring Function

```python
# In micromol_scorer.py, after VaccineScorer class:

class CustomScorer(VaccineScorer):
    """Domain-specific scoring."""
    
    def score_molecule(self, molecule_id, molecule_name, smiles, molecular_weight=None):
        # Call parent scorer
        base_score = super().score_molecule(molecule_id, molecule_name, smiles, molecular_weight)
        
        # Custom adjustments for your use case
        # E.g., if target is "cancer immunotherapy adjuvant":
        if self._is_cancer_related(molecule_name):
            base_score.overall_vaccine_score *= 1.2  # Boost score
            base_score.recommendation = "high" if base_score.overall_vaccine_score > 50 else "medium"
        
        return base_score
    
    def _is_cancer_related(self, name: str) -> bool:
        cancer_keywords = ["tumor", "cancer", "necrosis", "checkpoint"]
        return any(kw in name.lower() for kw in cancer_keywords)

# Usage:
scorer = CustomScorer()
score = scorer.score_molecule("my_mol", "Tumor necrosis factor", "CC(=O)O", 300)
```

### Export Network Visualization

```python
# Requires: pip install networkx plotly

import networkx as nx
import plotly.graph_objects as go

def visualize_metabolic_network(final_concentrations, reactions):
    """Create interactive network graph."""
    
    G = nx.DiGraph()
    
    # Add metabolite nodes
    for met_id, conc in final_concentrations.items():
        node_size = 20 + conc * 10  # Size by concentration
        G.add_node(met_id, size=node_size, type="metabolite")
    
    # Add reaction edges
    for rxn in reactions:
        for substrate in rxn.substrate_ids:
            for product in rxn.product_ids:
                G.add_edge(substrate, product, weight=1)
    
    # Layout
    pos = nx.spring_layout(G, k=2, iterations=50)
    
    # Plot
    edge_x = []
    edge_y = []
    for edge in G.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_x.append(x0)
        edge_x.append(x1)
        edge_x.append(None)
        edge_y.append(y0)
        edge_y.append(y1)
        edge_y.append(None)
    
    fig = go.Figure(
        data=[
            go.Scatter(
                x=edge_x, y=edge_y,
                mode='lines',
                line=dict(width=0.5, color='#888'),
                hoverinfo='none'
            ),
            go.Scatter(
                x=[pos[node][0] for node in G.nodes()],
                y=[pos[node][1] for node in G.nodes()],
                mode='markers+text',
                text=list(G.nodes()),
                textposition="top center",
                marker=dict(
                    size=[G.nodes[node].get('size', 10) for node in G.nodes()],
                    color='lightblue'
                )
            )
        ]
    )
    
    fig.show()
```

---

## TROUBLESHOOTING QUICK FIXES

| Problem | Solution |
|---------|----------|
| "No novel metabolites found" | Extend reactions (current models are simplified). Add more reactions from BiGG. |
| "Very low vaccine scores" | Adjust scoring thresholds in `VaccineScorer`. Update reference adjuvant database. |
| "Slow execution" | Use fewer workers, reduce duration, skip monocultures with `--skip-mono`. |
| "Memory error" | Use `--workers 2`, reduce num_steps in simulate(). |
| "RDKit errors" | Install: `pip install rdkit-pypi`. If it fails, simulator still works with fallback string similarity. |
| "Import error: micromol_core" | Ensure all 3 .py files in same directory. Check working directory: `pwd`. |

---

## NEXT IMMEDIATE ACTIONS

1. **Run the default simulation** (takes ~1 min):
   ```bash
   python run_simulator.py
   ```

2. **Check outputs**:
   ```bash
   cat micromol_results/micromol_summary.txt
   ```

3. **Customize for your use case**:
   - Add microbes (edit `create_default_microbes()`)
   - Update adjuvant database (edit `KnownMoleculeDatabase`)
   - Change scoring weights (edit `VaccineScorer.score_molecule()`)

4. **Run screen with your parameters**:
   ```bash
   python run_simulator.py --microbes YOUR_MICROBES --duration 48
   ```

5. **Read TECHNICAL_ROADMAP.md** for advanced features.

---

## PUBLICATIONS TO READ FIRST

**Essential Background:**
1. Orth et al. 2010 — "What is FBA?" ([Nature Biotechnology](https://doi.org/10.1038/nbt.1614))
2. Bordbar et al. 2015 — "Constraint-based models" ([Nature Rev Mol Cell Biol](https://doi.org/10.1038/nrm.2015.5))
3. Monk et al. 2014 — "iMM904: metabolic model" ([BMC Syst Biol](https://doi.org/10.1186/1752-0509-8-43))

**Vaccine Immunology:**
1. Rappuoli et al. 2016 — "Vaccines in the era of personalized medicine" ([Nature Reviews Immunology](https://doi.org/10.1038/nri.2016.42))
2. Pulendran & Ahmed 2011 — "Immunological mechanisms..." ([Cell](https://doi.org/10.1016/j.cell.2011.07.006))

**Microbiome + Immunity:**
1. Lynch & Pedersen 2016 — "The human microbiome..." ([NEJM](https://doi.org/10.1056/NEJMra1600266))

---

**You're ready to go!** 🚀🧬
