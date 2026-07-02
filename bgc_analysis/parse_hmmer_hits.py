#!/usr/bin/env python3
"""Parse an HMMER --tblout file, extract matching protein sequences from a FASTA,
and write candidate proteins and a CSV summary.
Usage:
  python parse_hmmer_hits.py --hits hits.tbl --proteins proteins.faa --out candidates
"""
import argparse
from pathlib import Path
import csv


def read_tblout(path):
    hits = []
    for line in open(path, 'r'):
        if line.startswith('#'):
            continue
        parts = line.strip().split()
        if len(parts) < 6:
            continue
        target_name = parts[0]
        query_name = parts[3]
        # store (target, query)
        hits.append((target_name, query_name))
    return hits


def read_fasta(path):
    seqs = {}
    cur = None
    for line in open(path, 'r'):
        line = line.rstrip('\n')
        if line.startswith('>'):
            cur = line[1:].split()[0]
            seqs[cur] = []
        else:
            if cur is not None:
                seqs[cur].append(line)
    for k in list(seqs.keys()):
        seqs[k] = '\n'.join(seqs[k])
    return seqs


def write_fasta(path, seqs):
    with open(path, 'w') as f:
        for sid, seq in seqs.items():
            f.write(f'>{sid}\n')
            # wrap at 80
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--hits', required=True)
    p.add_argument('--proteins', required=True)
    p.add_argument('--out', default='candidates')
    args = p.parse_args()

    hits = read_tblout(args.hits)
    seqs = read_fasta(args.proteins)

    # select unique target ids
    selected = {}
    for tid, qname in hits:
        if tid in seqs:
            selected[tid] = selected.get(tid, set())
            selected[tid].add(qname)

    out_prefix = Path(args.out)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    fasta_out = str(out_prefix) + '.faa'
    csv_out = str(out_prefix) + '.csv'

    write_fasta(fasta_out, {k: seqs[k] for k in selected.keys()})

    with open(csv_out, 'w', newline='') as csvf:
        writer = csv.writer(csvf)
        writer.writerow(['seq_id','hmm_queries'])
        for k, qset in selected.items():
            writer.writerow([k, ';'.join(sorted(qset))])

    print(f'Wrote {fasta_out} and {csv_out} ({len(selected)} candidates)')
