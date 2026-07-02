"""
Run media experiments: vary media formulations and run pairwise combinations, saving outputs per media.
"""
import argparse
import time
from pathlib import Path
from run_simulator import SimulationConfig, run_monocultures, run_all_combinations

MEDIA_PRESETS = {
    "glucose_only": {"glucose": 20.0},
    "glucose_plus_aa": {"glucose": 20.0, "alanine": 5.0, "glutamate": 5.0},
    "lactose": {"lactose": 20.0},
    "maltose": {"maltose": 20.0},
}


def run_for_media(media_name: str, duration: float = 24.0, serial: bool = True, workers: int = None):
    config = SimulationConfig()
    # Use expanded microbe list (all known)
    microbes = list(create_default_microbes().keys()) if 'create_default_microbes' in globals() else None


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--media', nargs='+', default=['glucose_only', 'glucose_plus_aa'])
    parser.add_argument('--duration', type=float, default=6.0)
    parser.add_argument('--serial', action='store_true')
    args = parser.parse_args()

    # Import here to avoid circular issues
    from micromol_core import create_default_microbes, create_metabolite_registry
    from run_simulator import SimulationConfig, run_monocultures, run_all_combinations

    base_out = Path('micromol_results')
    base_out.mkdir(exist_ok=True)

    for media in args.media:
        if media not in MEDIA_PRESETS:
            print(f"Unknown media {media}, skipping")
            continue
        out_dir = base_out / f"media_run_{media}"
        out_dir.mkdir(parents=True, exist_ok=True)

        print(f"Running media: {media} -> output: {out_dir}")
        # Prepare config
        config = SimulationConfig()
        config.microbe_ids = list(create_default_microbes().keys())
        config.simulation_duration = args.duration
        # Use pairwise combinations
        config.combination_sizes = [2]
        config.initial_substrate = None
        config.initial_concentration = None
        config.enable_parallel = not args.serial
        config.num_workers = None
        config.output_dir = out_dir

        # Set initial conditions from preset
        initial_conditions = MEDIA_PRESETS[media]

        # Run monocultures and combinations manually to control initial nutrients
        microbes = create_default_microbes()
        metabolites = create_metabolite_registry()
        from micromol_core import CombinationSimulator
        sim = CombinationSimulator(microbes, metabolites)

        # Monocultures
        print("  Running monocultures...")
        mono_results = []
        for m in config.microbe_ids:
            try:
                r = sim.simulate_monoculture(m, initial_conditions.copy(), duration=config.simulation_duration)
                r_entry = {"combination": m, "microbe_ids": [m], "novel_metabolites": r.get('novel_metabolites', []), "novel_metabolite_count": len(r.get('novel_metabolites', [])), "final_concentrations": r.get('final_concentrations', {}), "success": True}
                mono_results.append(r_entry)
            except Exception as e:
                mono_results.append({"combination": m, "error": str(e), "success": False})

        # Combinations
        print("  Running pairwise combinations...")
        # Reuse run_all_combinations but override initial conditions in worker by writing a small wrapper
        # Simpler: run serial here
        combos = []
        from itertools import combinations
        for combo in combinations(config.microbe_ids, 2):
            combos.append(list(combo))

        combo_results = []
        for combo in combos:
            try:
                r = sim.simulate_coculture(combo, initial_conditions.copy(), duration=config.simulation_duration)
                combo_results.append({"combination": " + ".join(combo), "microbe_ids": combo, "novel_metabolites": r.get('novel_metabolites', []), "novel_metabolite_count": len(r.get('novel_metabolites', [])), "scores": [], "final_concentrations": r.get('final_concentrations', {}), "success": True})
            except Exception as e:
                combo_results.append({"combination": " + ".join(combo), "error": str(e), "success": False})

        # Save results
        import json
        with open(out_dir / 'micromol_detailed.json', 'w') as f:
            json.dump(mono_results + combo_results, f, indent=2, default=str)

        # Write simple CSV summary
        import csv
        rows = []
        for r in combo_results:
            if r.get('success'):
                rows.append({
                    'Combination': r.get('combination'),
                    'Novel Metabolites Count': r.get('novel_metabolite_count', 0),
                    'Novel Metabolites': ", ".join(r.get('novel_metabolites', []))
                })
        if rows:
            with open(out_dir / 'micromol_results.csv', 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=rows[0].keys())
                writer.writeheader()
                writer.writerows(rows)

        print(f"  Finished media run: {media}")

    print("All media runs complete.")
