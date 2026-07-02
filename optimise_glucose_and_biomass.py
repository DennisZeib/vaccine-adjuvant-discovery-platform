"""
Optimise glucose concentration and biomass ratios for given microbe pairs.
Saves detailed JSON per run and an aggregate CSV summary.
"""
import json
import csv
from pathlib import Path
from micromol_core import CombinationSimulator, create_default_microbes, create_metabolite_registry
from micromol_scorer import create_default_scorer
import argparse


def parse_pair(s: str):
    # Accept formats like 'salmonella+clostridium' or 'a b'
    if '+' in s:
        return [p.strip() for p in s.split('+')]
    parts = s.split()
    if len(parts) == 2:
        return parts
    raise ValueError(f"Invalid pair format: {s}")


def run_sweep(pairs, concentrations, ratios, duration, outdir: Path):
    microbes = create_default_microbes()
    metabolites = create_metabolite_registry()
    simulator = CombinationSimulator(microbes, metabolites)
    scorer = create_default_scorer()

    outdir.mkdir(parents=True, exist_ok=True)
    summary_rows = []

    for pair in pairs:
        microA, microB = pair
        for conc in concentrations:
            initial = {"glucose": conc}
            # First sweep: equal biomass
            for r in ratios:
                # interpret r as ratio for microA relative to microB (microB fixed at 1.0)
                biomass = {microA: r, microB: 1.0}
                res = simulator.simulate_coculture(
                    microbe_ids=[microA, microB],
                    initial_nutrients=initial,
                    duration=duration,
                    biomass_ratios=biomass
                )
                # score top novel metabolite if exists
                novel = res.get("novel_metabolites", [])
                best_score = None
                best_met = None
                if novel:
                    for met in novel:
                        if met in metabolites:
                            s = scorer.score_molecule(met, metabolites[met].name, metabolites[met].smiles, metabolites[met].mw)
                            if best_score is None or s.overall_vaccine_score > best_score.overall_vaccine_score:
                                best_score = s
                                best_met = met
                final_conc = res.get("final_concentrations", {})
                lps_conc = final_conc.get("polysaccharide_lps", 0.0)

                # Save detailed JSON for this run
                fname = outdir / f"{microA}_{microB}_glc{int(conc)}_ratio{r}.json"
                with open(fname, 'w') as f:
                    json.dump({
                        "pair": [microA, microB],
                        "glucose": conc,
                        "ratio": biomass,
                        "result": res
                    }, f, indent=2)

                summary_rows.append({
                    "microA": microA,
                    "microB": microB,
                    "glucose_mM": conc,
                    "ratio_microA": r,
                    "lps_mM": lps_conc,
                    "best_novel_metabolite": best_met or "",
                    "best_score": best_score.overall_vaccine_score if best_score else ""
                })

    # Write summary CSV
    csvf = outdir / "optimisation_summary.csv"
    keys = ["microA", "microB", "glucose_mM", "ratio_microA", "lps_mM", "best_novel_metabolite", "best_score"]
    with open(csvf, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(summary_rows)

    print(f"Sweep completed. Results in {outdir}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Sweep glucose concentrations and biomass ratios for microbe pairs")
    parser.add_argument('--pair', action='append', required=True, help="Microbe pair, e.g. salmonella+clostridium (can specify multiple --pair)")
    parser.add_argument('--concentrations', default='10,20,50,100', help="Comma-separated glucose concentrations in mM")
    parser.add_argument('--ratios', default='0.1,0.2,0.5,1,2,5,10', help="Comma-separated biomass ratios (microA relative to microB)")
    parser.add_argument('--duration', type=float, default=48.0, help="Simulation duration in hours")
    parser.add_argument('--outdir', default='micromol_results/optimisation/glucose_biomass_sweep', help="Output directory")

    args = parser.parse_args()
    pairs = [parse_pair(p) for p in args.pair]
    concentrations = [float(x) for x in args.concentrations.split(',') if x.strip()]
    ratios = [float(x) for x in args.ratios.split(',') if x.strip()]
    outdir = Path(args.outdir)

    run_sweep(pairs, concentrations, ratios, args.duration, outdir)
