# MicroMol Combination Simulator - Setup & Usage Guide

## OVERVIEW

MicroMol is a **production-grade metabolic network simulator** for discovering vaccine candidates from microbe-molecule combinations.

**What it does:**
- Simulates metabolic networks using Monod kinetics
- Runs pairwise (and n-ary) microbe combinations
- Detects novel molecules produced only in co-cultures
- Scores molecules for vaccine potential (immunogenicity, safety, novelty)
- Parallelizes across CPU cores for 10,000+ combinations
- Outputs CSV, JSON, and human-readable reports

---

## INSTALLATION

### Prerequisites
- **Python 3.8+** (tested on 3.9, 3.10, 3.11)
- **pip** (Python package manager)
- **git** (optional, for cloning)

### Step 1: Clone/Extract Code

If you have a ZIP file:
```bash
unzip micromol_sim.zip
cd micromol_sim
```

Or create a folder manually:
```bash
mkdir micromol_sim
cd micromol_sim
# Copy the three .py files here
```

### Step 2: Create Virtual Environment (Recommended)

**On macOS/Linux:**
```bash
python3 -m venv venv
source venv/bin/activate
```

**On Windows (PowerShell):**
```powershell
python -m venv venv
.\venv\Scripts\Activate.ps1
```

**On Windows (Command Prompt):**
```cmd
python -m venv venv
venv\Scripts\activate.bat
```

### Step 3: Install Dependencies

Core dependencies (required):
```bash
pip install numpy pandas
```

**Optional but recommended** (for advanced features):
```bash
pip install rdkit-pypi
```

If `rdkit-pypi` fails, try:
```bash
pip install rdkit
```

**Full setup** (all features):
```bash
pip install numpy pandas rdkit-pypi
```

Verify installation:
```bash
python -c "import numpy, pandas; print('Core OK')"
python -c "from rdkit import Chem; print('RDKit OK')" # Optional
```

---

## QUICK START (5 Minutes)

### Run Default Simulation

In your `micromol_sim` directory:

```bash
python run_simulator.py
```

**Expected output:**
```
======================================================================
MICROMOL COMBINATION SIMULATOR
======================================================================
Microbes: ecoli, yeast, flu_virus
Duration: 24.0 hours
Initial substrate: glucose @ 10.0 mM
Parallel mode: True
Output directory: micromol_results
======================================================================

Running monoculture simulations...
  ecoli: 2 novel products
  yeast: 1 novel products
  flu_virus: 0 novel products

Generating combinations (max size 2)...
Total combinations to simulate: 3
Running 3 simulations with 8 workers...
Successful: 3, Failed: 0

Total time: 12.3s

Results saved to micromol_results/
  - micromol_results.csv (summary table)
  - micromol_detailed.json (full details)
  - micromol_summary.txt (human-readable report)
```

### Check Results

```bash
cat micromol_results/micromol_summary.txt
```

Or open in Excel/Google Sheets:
```bash
cat micromol_results/micromol_results.csv
```

---

## ADVANCED USAGE

### 1. Test with Different Microbes

Run only E. coli and yeast (no virus):
```bash
python run_simulator.py --microbes ecoli yeast
```

Add your own microbes (advanced):
```bash
python run_simulator.py --microbes ecoli yeast salmonella
```
(You'll need to add `salmonella` to `create_default_microbes()` first)

### 2. Longer Simulations

Simulate 48 hours:
```bash
python run_simulator.py --duration 48
```

### 3. Different Initial Substrate

Use acetate instead of glucose:
```bash
python run_simulator.py --substrate acetate --concentration 5.0
```

### 4. Serial Mode (Debugging)

Run without parallelization (useful for debugging):
```bash
python run_simulator.py --serial
```

### 5. Control Parallel Workers

Use 4 workers instead of auto-detect:
```bash
python run_simulator.py --workers 4
```

Use all CPUs - 1 (leave one core free):
```bash
python run_simulator.py --workers 6  # On 8-core system
```

### 6. Custom Output Directory

Save results elsewhere:
```bash
python run_simulator.py --output-dir ~/my_results/experiment_001
```

### 7. Skip Monocultures

Run only combinations (faster):
```bash
python run_simulator.py --skip-mono
```

### Complete Example Command

```bash
python run_simulator.py \
  --microbes ecoli yeast \
  --duration 36 \
  --substrate glucose \
  --concentration 15.0 \
  --workers 4 \
  --output-dir results/exp_20240123
```

---

## OUTPUT FILES

### `micromol_results.csv`
Spreadsheet format. Columns:
- **Combination**: Microbes in the co-culture
- **Novel Metabolites Count**: Number of unique products
- **Novel Metabolites**: List of metabolite IDs
- **Best Vaccine Score**: Highest vaccine potential (0-100)
- **Best Recommendation**: high/medium/low/excluded

### `micromol_detailed.json`
Full JSON with all scoring details:
- Final metabolite concentrations
- All novel molecules
- Complete scoring breakdown per molecule
- Immunogenicity, novelty, toxicity scores
- Rationale for each score

### `micromol_summary.txt`
Human-readable report:
- Monoculture results
- Top 5 candidate combinations
- Best vaccine candidates per combination

---

## CUSTOMIZATION

### Add a New Microbe

Edit `micromol_core.py` and add to `create_default_microbes()`:

```python
"mycobacterium": Microbe(
    id="mycobacterium",
    name="Mycobacterium tuberculosis",
    domain="bacteria",
    reactions=create_mycobacterium_reactions(),  # Define your own
    native_metabolites={"mycolic_acid"},
    essential_nutrients={"glucose"},
    max_growth_rate=0.2  # Slow grower
)
```

Define reactions:
```python
def create_mycobacterium_reactions() -> Dict[str, Reaction]:
    return {
        "m1_mycolic_acid_synthesis": Reaction(
            id="m1_mycolic_acid_synthesis",
            name="Acetyl-CoA → Mycolic acid",
            enzyme="Mycolic acid synthase",
            substrate_ids=["acetate"],
            product_ids=["mycolic_acid"],
            km_values={"acetate": 0.3},
            vmax=50.0
        ),
    }
```

Then run:
```bash
python run_simulator.py --microbes ecoli mycobacterium
```

### Add New Metabolites

Edit `create_metabolite_registry()` in `micromol_core.py`:

```python
"your_molecule": Metabolite(
    id="your_molecule",
    name="Human-readable name",
    smiles="CC(=O)O",  # SMILES string
    formula="C2H4O2",
    mw=60.05,
    source="nutrient"
),
```

Get SMILES from:
- PubChem (https://pubchem.ncbi.nlm.nih.gov/)
- ChemSpider (https://www.chemspider.com/)
- KEGG (https://www.kegg.jp/)

### Add Known Adjuvants

Edit `KnownMoleculeDatabase` in `micromol_scorer.py`:

```python
VACCINE_ADJUVANTS = {
    "my_adjuvant": {
        "smiles": "your_smiles_here",
        "name": "My Adjuvant",
        "class": "TLR agonist"
    },
    ...
}
```

---

## TROUBLESHOOTING

### Error: "No module named 'micromol_core'"

**Solution:** Make sure all three `.py` files are in the same directory:
```bash
ls -la *.py
# Should show: micromol_core.py, micromol_scorer.py, run_simulator.py
```

### Error: "RDKit not available"

**Solution:** Install RDKit (it's optional but recommended):
```bash
pip install rdkit-pypi
```

If that fails:
```bash
pip install rdkit
```

The simulator will work without RDKit but with reduced scoring accuracy.

### Memory Error on Large Simulations

**Solution:** Use fewer workers:
```bash
python run_simulator.py --workers 2
```

Or increase your system RAM / run on a machine with more memory.

### Slow Performance

**Solution:** 
- Use parallel mode (default): `python run_simulator.py`
- Run monocultures separately: `python run_simulator.py --skip-mono`
- Reduce duration: `python run_simulator.py --duration 12`

### Simulation Produces No Novel Metabolites

This is **normal** with the default simplified models. Real microbe genomes have thousands of reactions. To get more interesting results:
- Extend reactions with real BiGG database (next section)
- Add more microbes
- Use longer simulation times

---

## NEXT STEPS: Advanced Upgrades

### Level 1: Real Metabolic Models (BiGG Database)

Replace simplified reactions with genuine genome-scale metabolic models:

```python
pip install cobra

from cobra.io import read_sbml_model

# Load BiGG models
ecoli_model = read_sbml_model("e_coli_core.xml")  # Download from BiGG
```

Then integrate into `CombinationSimulator` via FBA (Flux Balance Analysis):
```python
def simulate_with_fba(self, microbe_ids, substrate):
    # Use COBRA's optimize() instead of Monod kinetics
    pass
```

Download free models: https://bigg.ucsd.edu/

### Level 2: Machine Learning Integration

Train a model to predict vaccine potential from sequence:

```python
from sklearn.ensemble import RandomForestClassifier

# Train on known vaccine adjuvants
X = [molecule_fingerprint(smiles) for smiles in training_data]
y = [is_adjuvant for _, is_adjuvant in training_data]
rf = RandomForestClassifier().fit(X, y)

# Predict on novel molecules
scores = rf.predict_proba(new_molecule_fingerprints)
```

### Level 3: Parallelization on Clusters

Scale to 100,000+ combinations using SLURM:

```bash
# Create sbatch script
cat > run_micromol.sbatch << EOF
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00

python run_simulator.py --workers 32 --duration 72
EOF

sbatch run_micromol.sbatch
```

### Level 4: Real Genomic Data Integration

Connect to reference genomes:

```python
from biopython import SeqIO

# Load actual genome
genome = SeqIO.read("escherichia_coli.fasta", "fasta")

# Predict genes → reactions → metabolites
genes = predict_genes(genome)  # Use Prodigal, GlimmerHMM, etc.
reactions = genes_to_reactions(genes)  # Map to BiGG
```

### Level 5: Immunogenicity Prediction (Deep Learning)

Replace heuristics with trained neural network:

```python
import torch
from transformers import AutoTokenizer, AutoModel

# Use pre-trained protein language model
tokenizer = AutoTokenizer.from_pretrained("facebook/esm2_t6_8M_UR50D")
model = AutoModel.from_pretrained("facebook/esm2_t6_8M_UR50D")

# Score peptides for MHC binding, T-cell epitopes, etc.
```

---

## EXAMPLE WORKFLOW: DISCOVERY SCENARIO

Suppose you're screening for **novel COVID-19 vaccine adjuvants**:

```bash
# 1. Create directory
mkdir covid_vaccine_screen
cd covid_vaccine_screen

# 2. Add SARS-CoV-2 to microbes (edit micromol_core.py)
# 3. Add known COVID antigens/adjuvants to KnownMoleculeDatabase (micromol_scorer.py)

# 4. Run screening: COVID virus + human immune cells
python run_simulator.py \
  --microbes sars_cov2 human_dendritic_cell macrophage \
  --duration 72 \
  --substrate glucose \
  --concentration 20.0 \
  --workers 16 \
  --output-dir covid_results

# 5. Check results
cat covid_results/micromol_summary.txt
# Look for: high vaccine scores, low toxicity, novel molecules

# 6. Validate top candidates (external wet-lab work)
# Test in animal models, clinical trials, etc.
```

---

## PUBLICATION / CITATION

If you use MicroMol in research, cite as:

```
MicroMol Combination Simulator v1.0
Real-time metabolic network analysis for vaccine discovery.
https://github.com/[your-repo]
```

---

## TECHNICAL DETAILS

### Simulation Engine

**Kinetic Model:** Monod kinetics with Michaelis-Menten saturation
```
v = Vmax × [S] / (Km + [S])
```

**Time Integration:** Explicit Euler method (dt = 0.1h)

**Metabolite Exchange:** Passive diffusion (secretion) between microbes

### Scoring Algorithm

**Vaccine Potential Score (0-100):**
```
score = (immunogenicity + novelty + 0.8×adjuvant_sim) 
        × (1 - toxicity_risk) 
        × physicochemical_score 
        × 100/3
```

**Immunogenicity Heuristic:**
- Molecular weight (200-1500 Da optimal)
- Aromatic rings, heteroatoms
- Polar surface area (20-150 Ų)
- Charge, hydrogen bonding

**Novelty Score:**
- 1 - (average similarity to all known molecules)
- Similarity via Morgan fingerprints (RDKit)
- Tanimoto coefficient

**Toxicity Risk:**
- Similarity to known toxins
- Lipinski violations (MW > 500, LogP > 5, etc.)
- Known dangerous functional groups

---

## FAQ

**Q: Can I run this on Windows?**
A: Yes. Install Python 3.8+ from python.org. Then follow the Windows instructions in "Installation."

**Q: Can I add my own microbes?**
A: Yes. Edit `create_default_microbes()` in `micromol_core.py` and add a new `Microbe` object.

**Q: How do I integrate real genomic data?**
A: Download SBML models from BiGG (https://bigg.ucsd.edu/), then use COBRA library for FBA instead of Monod kinetics.

**Q: What's the maximum number of combinations?**
A: Practically unlimited. With 10 microbes and triplet combinations, that's ~120 simulations. Modern hardware can handle 100,000+ with parallelization.

**Q: Can this predict clinical vaccine efficacy?**
A: No. This is a computational screening tool. Validated candidates must be tested in vitro and in vivo.

---

## SUPPORT

- Check `micromol_results/micromol_summary.txt` for detailed logs
- Enable debug logging:
  ```python
  logging.basicConfig(level=logging.DEBUG)
  ```
- File an issue on GitHub with:
  - Python version
  - Error message
  - Command used
  - System specs (CPU, RAM)

---

**Happy discovering!** 🧬🔬
