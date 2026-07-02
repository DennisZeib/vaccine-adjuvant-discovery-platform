#!/usr/bin/env python3
"""Generate candidate proteins from a proteins FASTA using simple heuristics.
Writes <out>.faa and <out>.csv compatible with `score_candidates.py`.

Usage:
  python generate_candidates_from_proteins.py --proteins genome/proteins.faa --out genome/candidates --minlen 300 --top 200
"""
import argparse
from pathlib import Path
import heapq


def read_fasta(path):
    seqs = {}
    cur = None
    with open(path, 'r') as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                cur = line[1:].split()[0]
                seqs[cur] = []
            else:
                if cur is not None:
                    seqs[cur].append(line)
    for k in list(seqs.keys()):
        seqs[k] = ''.join(seqs[k])
    return seqs


def write_fasta(path, seqs):
    with open(path, 'w') as f:
        for sid, seq in seqs.items():
            f.write(f'>{sid}\n')
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--proteins', required=True)
    p.add_argument('--out', default='genome/candidates')
    p.add_argument('--minlen', type=int, default=300)
    p.add_argument('--top', type=int, default=200, help='Limit to top N longest proteins')
    args = p.parse_args()

    proteins = Path(args.proteins)
    if not proteins.exists():
        print(f'Proteins FASTA not found: {proteins}')
        return

    seqs = read_fasta(proteins)
    # select by length
    heap = []
    for sid, seq in seqs.items():
        ln = len(seq)
        if ln >= args.minlen:
            import heapq
            heapq.heappush(heap, (-ln, sid))

    selected = []
    while heap and len(selected) < args.top:
        _, sid = heapq.heappop(heap)
        selected.append(sid)

    out_pref = Path(args.out)
    out_pref.parent.mkdir(parents=True, exist_ok=True)
    fasta_out = str(out_pref) + '.faa'
    csv_out = str(out_pref) + '.csv'

    write_fasta(fasta_out, {k: seqs[k] for k in selected})

    with open(csv_out, 'w') as f:
        f.write('seq_id,hmm_queries\n')
        for sid in selected:
            f.write(f'{sid},LENGTH>={args.minlen}\n')

    print(f'Wrote {fasta_out} and {csv_out} ({len(selected)} candidates)')


if __name__ == '__main__':
    main()
