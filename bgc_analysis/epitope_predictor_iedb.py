#!/usr/bin/env python3
"""IEDB epitope predictor helper

Reads `bgc_analysis/reports/top_candidates.csv` (top N proteins), extracts sequences
from `bgc_analysis/genome/candidates.faa`, queries the IEDB MHC I API, and writes
`bgc_analysis/reports/top10_epitopes.csv` with predicted binders.

Usage:
  python epitope_predictor_iedb.py --top 10 --alleles HLA-A02:01,HLA-B07:02 --lengths 9 --dry-run

Dry-run performs sequence extraction and writes a placeholder output without network calls.
"""
import argparse
import csv
import sys
from pathlib import Path

BASE = Path(__file__).resolve().parent
REPORTS = BASE / 'reports'
GENOME = BASE / 'genome'
TOP_CSV = REPORTS / 'top_candidates.csv'
CANDIDATES_FASTA = GENOME / 'candidates.faa'
OUT_CSV = REPORTS / 'top10_epitopes.csv'

IEDB_MHCI_URL = 'http://tools-cluster-interface.iedb.org/tools_api/mhci'


def read_top_ids(path, top=10):
    ids = []
    with open(path, newline='') as fh:
        r = csv.DictReader(fh)
        for i, row in enumerate(r):
            if i >= top:
                break
            ids.append(row['Protein_ID'])
    return ids


def read_fasta(path):
    seqs = {}
    with open(path) as fh:
        name = None
        parts = []
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            if line.startswith('>'):
                if name:
                    seqs[name] = ''.join(parts)
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
        if name:
            seqs[name] = ''.join(parts)
    return seqs


def call_iedb_mhci(sequence, alleles, lengths, method='netmhcpan'):
    """Call IEDB mhci API. Returns list of dicts with keys: allele, start, length, peptide, ic50, rank"""
    import requests
    data = {
        'sequence_text': sequence,
        'alleles': ','.join(alleles),
        'length': ','.join(str(l) for l in lengths),
        'method': method,
    }
    resp = requests.post(IEDB_MHCI_URL, data=data, timeout=120)
    resp.raise_for_status()
    # API returns TSV with header starting with '#'
    lines = [l for l in resp.text.splitlines() if not l.startswith('#') and l.strip()]
    results = []
    for line in lines:
        # expected TSV columns: allele\tstart\tend\tpeptide\thic50\trank\n (variable)
        cols = line.split('\t')
        try:
            allele = cols[0]
            start = int(cols[1])
            peptide = cols[3]
            ic50 = float(cols[4]) if cols[4] not in ('NA','-') else None
            rank = float(cols[5]) if len(cols) > 5 and cols[5] not in ('NA','-') else None
            results.append({'allele': allele, 'start': start, 'peptide': peptide, 'ic50': ic50, 'rank': rank})
        except Exception:
            continue
    return results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--top', type=int, default=10)
    parser.add_argument('--alleles', type=str, default='HLA-A02:01')
    parser.add_argument('--lengths', type=str, default='9')
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()

    if not TOP_CSV.exists():
        print(f"Missing top candidates CSV: {TOP_CSV}")
        sys.exit(1)
    if not CANDIDATES_FASTA.exists():
        print(f"Missing candidates FASTA: {CANDIDATES_FASTA}")
        sys.exit(1)

    ids = read_top_ids(TOP_CSV, top=args.top)
    seqs = read_fasta(CANDIDATES_FASTA)

    alleles = [a.strip() for a in args.alleles.split(',') if a.strip()]
    lengths = [int(x) for x in args.lengths.split(',')]

    out_rows = []

    for pid in ids:
        seq = seqs.get(pid)
        if seq is None:
            print(f"Warning: sequence {pid} not found in {CANDIDATES_FASTA}")
            continue
        print(f"Processing {pid}: length {len(seq)}")
        if args.dry_run:
            # create a fake hit for testing
            out_rows.append({'Protein_ID': pid, 'Allele': alleles[0], 'Peptide': seq[:9], 'Start': 1, 'Length': 9, 'IC50': 'NA', 'Percentile_Rank': 'NA'})
            continue
        try:
            hits = call_iedb_mhci(seq, alleles, lengths)
        except Exception as e:
            print(f"IEDB API call failed for {pid}: {e}")
            continue
        for h in hits:
            out_rows.append({'Protein_ID': pid, 'Allele': h['allele'], 'Peptide': h['peptide'], 'Start': h['start'], 'Length': len(h['peptide']), 'IC50': h['ic50'] if h['ic50'] is not None else 'NA', 'Percentile_Rank': h['rank'] if h['rank'] is not None else 'NA'})

    # write output CSV
    REPORTS.mkdir(parents=True, exist_ok=True)
    with open(OUT_CSV, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['Protein_ID', 'Allele', 'Peptide', 'Start', 'Length', 'IC50', 'Percentile_Rank'])
        w.writeheader()
        for r in out_rows:
            w.writerow(r)

    print(f"Wrote epitopes to {OUT_CSV}")
    if args.dry_run:
        print("Dry-run: no network calls were made. Re-run without --dry-run to call IEDB API.")


if __name__ == '__main__':
    main()
