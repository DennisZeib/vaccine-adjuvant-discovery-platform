"""
MicroMol Simulator - Main Entry Point
Orchestrates combination screening with parallel processing and comprehensive output.
"""

import json
import csv
import time
import logging
from typing import List, Dict, Tuple, Any
from dataclasses import asdict
from multiprocessing import Pool, cpu_count
from pathlib import Path
from itertools import combinations
import argparse
import sys

from micromol_core import (
    CombinationSimulator, 
    MetabolicNetwork,
    create_default_microbes,
    create_metabolite_registry,
    Reaction,
    Microbe
)
from micromol_scorer import create_default_scorer, VaccinePotentialScore

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s | %(levelname)s | %(message)s'
)
logger = logging.getLogger(__name__)


class SimulationConfig:
    """Configuration for simulation run."""
    
    def __init__(self):
        self.microbe_ids = ["ecoli", "yeast", "flu_virus"]
        self.combination_sizes = [2]  # Pairwise combinations only
        self.simulation_duration = 24.0  # hours
        self.initial_substrate = "glucose"
        self.initial_concentration = 10.0  # mM
        self.num_workers = None  # Auto-detect CPUs
        self.enable_parallel = True
        self.output_dir = Path("micromol_results")
        self.run_monocultures = True
        self.run_all_interactions = True
        self.biomass_ratios = None


def run_single_combination(args: Tuple) -> Dict[str, Any]:
    """
    Worker function: simulate a single combination.
    Must be at module level for multiprocessing.
    
    Args:
        args: Tuple of (microbe_ids, initial_conditions, duration)
    
    Returns:
        Result dict with simulation and scoring data
    """
    microbe_ids, initial_conditions, duration, biomass_ratios = args
    
    try:
        # Initialize fresh simulator for this worker
        microbes = create_default_microbes()
        metabolites = create_metabolite_registry()
        simulator = CombinationSimulator(microbes, metabolites)
        scorer = create_default_scorer()
        
        # Run simulation
        combo_name = " + ".join(microbe_ids)
        result = simulator.simulate_coculture(
            microbe_ids=microbe_ids,
            initial_nutrients=initial_conditions,
            duration=duration,
            biomass_ratios=biomass_ratios
        )
        
        # Extract novel metabolites
        novel_mets = result.get("novel_metabolites", [])
        
        # Score novel molecules
        scores = []
        for met_id in novel_mets:
            if met_id in metabolites:
                met = metabolites[met_id]
                score = scorer.score_molecule(
                    molecule_id=met_id,
                    molecule_name=met.name,
                    smiles=met.smiles,
                    molecular_weight=met.mw
                )
                scores.append(asdict(score))
        
        return {
            "combination": combo_name,
            "microbe_ids": microbe_ids,
            "novel_metabolites": novel_mets,
            "novel_metabolite_count": len(novel_mets),
            "scores": scores,
            "final_concentrations": result.get("final_concentrations", {}),
            "success": True
        }
    
    except Exception as e:
        logger.error(f"Error simulating {microbe_ids}: {e}")
        return {
            "combination": " + ".join(microbe_ids),
            "microbe_ids": microbe_ids,
            "novel_metabolites": [],
            "novel_metabolite_count": 0,
            "scores": [],
            "error": str(e),
            "success": False
        }


def run_monocultures(config: SimulationConfig) -> List[Dict[str, Any]]:
    """Run simulations of single microbes."""
    logger.info("Running monoculture simulations...")
    
    results = []
    microbes = create_default_microbes()
    metabolites = create_metabolite_registry()
    simulator = CombinationSimulator(microbes, metabolites)
    
    for microbe_id in config.microbe_ids:
        try:
            initial = {config.initial_substrate: config.initial_concentration}
            result = simulator.simulate_monoculture(
                microbe_id=microbe_id,
                initial_nutrients=initial,
                duration=config.simulation_duration,
                biomass_ratios=config.biomass_ratios
            )
            
            result_entry = {
                "combination": microbe_id,
                "microbe_ids": [microbe_id],
                "novel_metabolites": result.get("novel_metabolites", []),
                "novel_metabolite_count": len(result.get("novel_metabolites", [])),
                "final_concentrations": result.get("final_concentrations", {}),
                "type": "monoculture",
                "success": True
            }
            results.append(result_entry)
            logger.info(f"  {microbe_id}: {result_entry['novel_metabolite_count']} novel products")
        
        except Exception as e:
            logger.error(f"Error in monoculture {microbe_id}: {e}")
            results.append({
                "combination": microbe_id,
                "microbe_ids": [microbe_id],
                "error": str(e),
                "success": False
            })
    
    return results


def generate_combinations(microbe_ids: List[str], 
                         max_size: int = 2) -> List[List[str]]:
    """Generate all combinations of microbe IDs up to max_size."""
    combos = []
    for size in range(2, min(max_size + 1, len(microbe_ids) + 1)):
        combos.extend(combinations(microbe_ids, size))
    return [list(c) for c in combos]


def run_all_combinations(config: SimulationConfig) -> List[Dict[str, Any]]:
    """Run all microbe combination simulations (with parallelization)."""
    logger.info(f"Generating combinations (max size {max(config.combination_sizes)})...")
    
    # Generate all combinations
    all_combos = []
    for size in config.combination_sizes:
        all_combos.extend(generate_combinations(config.microbe_ids, size))
    
    logger.info(f"Total combinations to simulate: {len(all_combos)}")
    
    # Prepare worker arguments
    initial_conditions = {config.initial_substrate: config.initial_concentration}
    worker_args = [
        (combo, initial_conditions, config.simulation_duration, config.biomass_ratios)
        for combo in all_combos
    ]
    
    # Run in parallel or serial
    results = []
    if config.enable_parallel:
        num_workers = config.num_workers or cpu_count()
        logger.info(f"Running {len(all_combos)} simulations with {num_workers} workers...")
        
        with Pool(num_workers) as pool:
            results = pool.map(run_single_combination, worker_args)
    else:
        logger.info(f"Running {len(all_combos)} simulations (serial mode)...")
        results = [run_single_combination(args) for args in worker_args]
    
    # Filter out errors
    successful = [r for r in results if r.get("success", False)]
    failed = [r for r in results if not r.get("success", False)]
    
    logger.info(f"Successful: {len(successful)}, Failed: {len(failed)}")
    
    return successful, failed


def write_results_csv(results: List[Dict[str, Any]], 
                     filename: str = "micromol_results.csv"):
    """Write combination results to CSV."""
    logger.info(f"Writing results to {filename}...")
    
    rows = []
    for result in results:
        if not result.get("success", False):
            continue
        
        combination = result.get("combination", "")
        novel_count = result.get("novel_metabolite_count", 0)
        novel_mets = ", ".join(result.get("novel_metabolites", []))
        
        # Get best vaccine score if available
        best_score = None
        best_recommendation = "N/A"
        if result.get("scores"):
            best = max(result["scores"], key=lambda s: s.get("overall_vaccine_score", 0))
            best_score = best.get("overall_vaccine_score", 0)
            best_recommendation = best.get("recommendation", "N/A")
        
        rows.append({
            "Combination": combination,
            "Novel Metabolites Count": novel_count,
            "Novel Metabolites": novel_mets,
            "Best Vaccine Score": best_score if best_score is not None else "N/A",
            "Best Recommendation": best_recommendation
        })
    
    # Sort by vaccine score (descending)
    rows.sort(
        key=lambda r: r["Best Vaccine Score"] 
        if isinstance(r["Best Vaccine Score"], (int, float)) else 0,
        reverse=True
    )
    
    # Write CSV
    if rows:
        with open(filename, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)
        logger.info(f"  Wrote {len(rows)} rows to {filename}")
    else:
        logger.warning(f"  No successful results to write")


def write_detailed_json(results: List[Dict[str, Any]], 
                       filename: str = "micromol_detailed.json"):
    """Write detailed results (with all scores) to JSON."""
    logger.info(f"Writing detailed results to {filename}...")
    
    # Clean results for JSON serialization
    clean_results = []
    for r in results:
        if not r.get("success", False):
            continue
        clean_r = dict(r)
        # Keep scores as-is (already dicts)
        clean_results.append(clean_r)
    
    with open(filename, "w") as f:
        json.dump(clean_results, f, indent=2, default=str)
    
    logger.info(f"  Wrote {len(clean_results)} detailed results")


def write_summary_report(all_results: Dict[str, List], 
                        filename: str = "micromol_summary.txt"):
    """Write a human-readable summary report."""
    logger.info(f"Writing summary report to {filename}...")
    
    with open(filename, "w") as f:
        f.write("=" * 70 + "\n")
        f.write("MICROMOL COMBINATION SIMULATOR - SUMMARY REPORT\n")
        f.write("=" * 70 + "\n\n")
        
        f.write(f"Timestamp: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        # Monoculture summary
        mono_results = all_results.get("monoculture", [])
        f.write(f"MONOCULTURE RESULTS ({len(mono_results)} runs):\n")
        f.write("-" * 70 + "\n")
        for result in mono_results:
            f.write(f"  {result.get('combination', 'N/A')}\n")
            f.write(f"    Novel metabolites: {result.get('novel_metabolite_count', 0)}\n")
        
        f.write("\n")
        
        # Combination summary
        combo_results = all_results.get("combinations", [])[0]  # successful
        f.write(f"COMBINATION RESULTS ({len(combo_results)} combinations):\n")
        f.write("-" * 70 + "\n")
        
        # Top candidates
        top_candidates = sorted(
            combo_results,
            key=lambda r: max(
                (s.get("overall_vaccine_score", 0) for s in r.get("scores", [])),
                default=0
            ),
            reverse=True
        )[:5]
        
        f.write("TOP 5 CANDIDATES (by vaccine potential):\n")
        for i, result in enumerate(top_candidates, 1):
            f.write(f"\n  {i}. {result.get('combination', 'N/A')}\n")
            f.write(f"     Novel metabolites: {result.get('novel_metabolite_count', 0)}\n")
            if result.get("scores"):
                best = max(result["scores"], key=lambda s: s.get("overall_vaccine_score", 0))
                f.write(f"     Best candidate: {best.get('molecule_name', 'N/A')}\n")
                f.write(f"     Vaccine score: {best.get('overall_vaccine_score', 0):.1f}/100\n")
                f.write(f"     Recommendation: {best.get('recommendation', 'N/A')}\n")
        
        f.write("\n" + "=" * 70 + "\n")
        f.write("END OF REPORT\n")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="MicroMol Combination Simulator - Metabolic network analysis and vaccine discovery"
    )
    parser.add_argument(
        "--microbes",
        nargs="+",
        default=["ecoli", "yeast", "flu_virus"],
        help="Microbe IDs to use (default: ecoli yeast flu_virus)"
    )
    parser.add_argument(
        "--duration",
        type=float,
        default=24.0,
        help="Simulation duration in hours (default: 24)"
    )
    parser.add_argument(
        "--substrate",
        default="glucose",
        help="Initial substrate (default: glucose)"
    )
    parser.add_argument(
        "--concentration",
        type=float,
        default=10.0,
        help="Initial substrate concentration in mM (default: 10)"
    )
    parser.add_argument(
        "--serial",
        action="store_true",
        help="Run in serial mode (no parallelization)"
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=None,
        help="Number of parallel workers (default: CPU count)"
    )
    parser.add_argument(
        "--output-dir",
        default="micromol_results",
        help="Output directory (default: micromol_results)"
    )
    parser.add_argument(
        "--skip-mono",
        action="store_true",
        help="Skip monoculture simulations"
    )
    parser.add_argument(
        "--biomass",
        default=None,
        help="Initial biomass ratios as comma-separated 'microbe=ratio' pairs (e.g. salmonella=1.0,lactobacillus=0.5)"
    )
    
    args = parser.parse_args()
    
    # Configure
    config = SimulationConfig()
    config.microbe_ids = args.microbes
    config.simulation_duration = args.duration
    config.initial_substrate = args.substrate
    config.initial_concentration = args.concentration
    config.enable_parallel = not args.serial
    config.num_workers = args.workers
    config.output_dir = Path(args.output_dir)
    config.run_monocultures = not args.skip_mono
    # Parse biomass ratios if provided
    def parse_biomass_arg(s: str) -> Dict[str, float]:
        out = {}
        if not s:
            return out
        parts = s.split(',')
        for p in parts:
            if '=' in p:
                k, v = p.split('=', 1)
            elif ':' in p:
                k, v = p.split(':', 1)
            else:
                continue
            try:
                out[k.strip()] = float(v)
            except ValueError:
                logger.warning(f"Invalid biomass value for {k}: {v}")
        return out

    config.biomass_ratios = parse_biomass_arg(args.biomass)
    
    # Create output directory
    config.output_dir.mkdir(parents=True, exist_ok=True)
    
    logger.info("=" * 70)
    logger.info("MICROMOL COMBINATION SIMULATOR")
    logger.info("=" * 70)
    logger.info(f"Microbes: {', '.join(config.microbe_ids)}")
    logger.info(f"Duration: {config.simulation_duration} hours")
    logger.info(f"Initial substrate: {config.initial_substrate} @ {config.initial_concentration} mM")
    logger.info(f"Parallel mode: {config.enable_parallel}")
    logger.info(f"Output directory: {config.output_dir}")
    logger.info("=" * 70 + "\n")
    
    # Run simulations
    start_time = time.time()
    all_results = {}
    
    # Monocultures
    if config.run_monocultures:
        mono_results = run_monocultures(config)
        all_results["monoculture"] = mono_results
    
    # Combinations
    combo_results, failed_results = run_all_combinations(config)
    all_results["combinations"] = [combo_results, failed_results]
    
    elapsed = time.time() - start_time
    logger.info(f"\nTotal time: {elapsed:.1f}s\n")
    
    # Write outputs
    output_dir = config.output_dir
    write_results_csv(combo_results, str(output_dir / "micromol_results.csv"))
    write_detailed_json(combo_results, str(output_dir / "micromol_detailed.json"))
    write_summary_report(all_results, str(output_dir / "micromol_summary.txt"))
    
    logger.info(f"\nResults saved to {output_dir}/")
    logger.info("  - micromol_results.csv (summary table)")
    logger.info("  - micromol_detailed.json (full details)")
    logger.info("  - micromol_summary.txt (human-readable report)")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
