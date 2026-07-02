#!/usr/bin/env python3
import pandas as pd
import sys
import os
from pathlib import Path

# allow importing micromol_scorer from repo root
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from micromol_scorer import create_default_scorer


def infer_domain_type(hmm_queries: str) -> str:
    q = (hmm_queries or "").upper()
    if 'PKS' in q or 'KETO' in q or 'KS_' in q:
        return 'PKS'
    if 'NRPS' in q or 'CONDENS' in q or 'A_DOMAIN' in q or 'AMP' in q:
        return 'NRPS'
    # fallback
    return 'PKS/NRPS'


def main():
    base = Path(__file__).parent
    cand_csv = base / 'genome' / 'candidates.csv'
    if not cand_csv.exists():
        print(f'Candidate CSV not found: {cand_csv}. Run parse_hmmer_hits.py first.')
        return

    df = pd.read_csv(cand_csv)
    # Support both expected column names
    possible_id_cols = ['seq_id', 'protein_id', 'seqid', 'id']
    id_col = next((c for c in df.columns if c in possible_id_cols), None)
    query_col = next((c for c in df.columns if 'hmm' in c.lower() or 'query' in c.lower() or 'hmm_queries'==c), None)
    if id_col is None:
        print('No sequence id column found in candidates.csv. Columns:', df.columns.tolist())
        return
    if query_col is None:
        # create empty
        df['hmm_queries'] = ''
        query_col = 'hmm_queries'

    scorer = create_default_scorer()
    results = []
    for _, row in df.iterrows():
        seq_id = row[id_col]
        hmm_q = row.get(query_col, '')
        domain = infer_domain_type(hmm_q)
        if domain == 'PKS':
            placeholder_smiles = 'CC(C)C1CC(=O)OC1'  # simple lactone
        elif domain == 'NRPS':
            placeholder_smiles = 'C1CC(=O)NC1'      # simple lactam
        else:
            placeholder_smiles = 'CC(C)C1CC(=O)OC1'  # default to PKS-like

        score = scorer.score_molecule(str(seq_id), f'Novel {domain} product', placeholder_smiles, 500)
        results.append({
            'Protein_ID': seq_id,
            'Domain_Type': domain,
            'Vaccine_Score': round(score.overall_vaccine_score, 3),
            'Recommendation': score.recommendation,
            'Rationale': '; '.join(score.rationale[:3])
        })

    out_df = pd.DataFrame(results).sort_values('Vaccine_Score', ascending=False)
    out_path = base / 'genome' / 'candidate_scores.csv'
    out_df.to_csv(out_path, index=False)
    print(f'Scores saved to {out_path}')
    print(out_df.head(10).to_string(index=False))


if __name__ == '__main__':
    main()
