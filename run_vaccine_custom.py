"""
Helper script: run a short targeted simulation using custom microbes that produce vaccine-relevant molecules.
"""
import sys
from pathlib import Path

# Ensure workspace path is on sys.path
root = Path(__file__).resolve().parent
sys.path.insert(0, str(root))

from run_simulator import main

if __name__ == "__main__":
    # Arguments: small 6-hour serial run with higher initial glucose
    sys.argv = [
        "run_simulator.py",
        "--microbes",
        "ecoli",
        "yeast",
        "bacillus",
        "salmonella",
        "lactobacillus",
        "--duration",
        "6",
        "--concentration",
        "20",
        "--serial",
        "--output-dir",
        "micromol_results/custom_run",
        "--skip-mono"
    ]
    sys.exit(main())
