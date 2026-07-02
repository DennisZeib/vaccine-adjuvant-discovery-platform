"""
Run expanded screen using all microbes defined in micromol_core.create_default_microbes()
"""
import sys
from pathlib import Path
root = Path(__file__).resolve().parent
sys.path.insert(0, str(root))

from micromol_core import create_default_microbes
from run_simulator import main

if __name__ == '__main__':
    microbes = list(create_default_microbes().keys())
    args = [
        'run_simulator.py',
        '--microbes',
    ] + microbes + [
        '--duration', '24',
        '--concentration', '20',
        '--workers', '4',
        '--output-dir', 'micromol_results/expanded_run'
    ]
    print('Running with microbes:', microbes)
    import run_simulator as rs
    # call main with args
    import sys
    sys.argv = args
    rs.main()
