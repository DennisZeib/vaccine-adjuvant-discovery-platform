#!/usr/bin/env python3
"""
discover_novel_bgc.py

Download a Streptomyces genome, attempt antiSMASH if installed, otherwise
create a placeholder PKS product SMILES and score it with micromol_scorer.

Usage: python discover_novel_bgc.py
"""
import os
import sys
from pathlib import Path
import shutil
import subprocess
import json
import urllib.request

BASE_DIR = Path(__file__).parent
GENOMES_DIR = BASE_DIR / "genomes"
GENOMES_DIR.mkdir(exist_ok=True)

# Example Streptomyces avermitilis genome (NCBI RefSeq)
URL = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/765/GCF_000009765.2_ASM976v2/GCF_000009765.2_ASM976v2_genomic.fna.gz"
OUT_GZ = GENOMES_DIR / "GCF_000009765.2_genomic.fna.gz"
OUT_FASTA = GENOMES_DIR / "GCF_000009765.2_genomic.fna"

def download_genome():
    if OUT_FASTA.exists():
        print(f"Genome already downloaded: {OUT_FASTA}")
        return
    print(f"Downloading genome to: {OUT_GZ}")
    urllib.request.urlretrieve(URL, OUT_GZ)
    try:
        import gzip
        with gzip.open(OUT_GZ, 'rb') as f_in, open(OUT_FASTA, 'wb') as f_out:
            f_out.write(f_in.read())
        print(f"Decompressed to: {OUT_FASTA}")
    except Exception as e:
        print("Failed to decompress automatically; you can unzip manually.")

def run_antismash(fasta_path: Path):
    antismash = shutil.which('antismash') or shutil.which('run_antismash')
    if not antismash:
        print("antismash not found in PATH; skipping local antiSMASH run.")
        return None
    outdir = GENOMES_DIR / "antismash_output"
    cmd = [antismash, str(fasta_path), '--output-dir', str(outdir)]
    print("Running antiSMASH (this may take several minutes)...")
    subprocess.run(cmd, check=False)
    if outdir.exists():
        print(f"antiSMASH output written to: {outdir}")
        return outdir
    return None

def fallback_score_placeholder():
    # Create a generic macrolactone-like placeholder SMILES for a polyketide
    placeholder_smiles = 'CCCC(=O)OCC1OC(=O)CCCCCC1'  # small macrolactone-like fragment
    try:
        from micromol_scorer import create_default_scorer
    except Exception as e:
        print("micromol_scorer not importable:", e)
        return None
    scorer = create_default_scorer()
    print("Scoring placeholder PKS product (novel_pks_product)...")
    score = scorer.score_molecule(
        molecule_id='novel_pks_product',
        molecule_name='Novel PKS product (placeholder)',
        smiles=placeholder_smiles,
        molecular_weight=500.0
    )
    out = {
        'molecule_id': 'novel_pks_product',
        'smiles': placeholder_smiles,
        'toxicity_risk': getattr(score, 'toxicity_risk', None),
        'immunogenicity_score': getattr(score, 'immunogenicity_score', None),
        'overall_vaccine_score': getattr(score, 'overall_vaccine_score', None)
    }
    summary_file = BASE_DIR / 'novel_candidate_summary.json'
    with open(summary_file, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"Saved candidate summary: {summary_file}")
    print(json.dumps(out, indent=2))
    return out

def main():
    download_genome()
    antismash_out = run_antismash(OUT_FASTA)
    if antismash_out:
        print("antiSMASH ran; you should inspect clusters for 'no known product' clusters in the output directory.")
        # For now, stop here; user can pick a cluster for structure prediction.
        return
    # fallback: generate and score placeholder
    res = fallback_score_placeholder()
    if res:
        print("Placeholder candidate scored. If you want, I can (A) attempt automatic structure generation from a selected BGC (requires antiSMASH/PRISM), or (B) run deeper searches across AGORA models to find enzymes that match PKS domains.")

if __name__ == '__main__':
    main()
