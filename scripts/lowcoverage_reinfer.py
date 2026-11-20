"""Infer from a previously generated synthetic cohort, retaining fixed truth."""
import argparse
import json
from pathlib import Path

import numpy as np

from ghist_lowcoverage.pipeline import run_pipeline, sha256, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--gl-model', type=int, choices=[1, 2], required=True)
    parser.add_argument('--angsd', type=Path, required=True)
    parser.add_argument('--realsfs', type=Path, required=True)
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--restarts', type=int, default=8)
    args = parser.parse_args()
    simulation = json.loads((args.case / 'synthetic/simulation.json').read_text())
    if simulation['final_data_access'] or simulation['testing_data_access']:
        raise ValueError('this harness is for independent synthetic data only')
    bam_dir = args.case / ('remapped' if (args.case / 'remapped').is_dir() else 'synthetic')
    result = run_pipeline(bam_dir, bam_dir / 'reference.fasta', simulation['task'], args.output,
                          args.angsd, args.realsfs, threads=args.threads, expected_length=simulation['length'],
                          restarts=args.restarts, genotype_likelihood_model=args.gl_model)
    errors = {k: result['targets'][k] / value - 1 for k, value in simulation['truth'].items()}
    write_json(args.output / 'scored.json', {'task': simulation['task'], 'gl_model': args.gl_model,
        'source_sha256': sha256(Path(__file__)), 'truth': simulation['truth'], 'estimate': result['targets'],
        'relative_error': errors, 'diagnostic_rrmse': float(np.sqrt(np.mean(np.square(list(errors.values()))))),
        'numerically_valid': result['numerically_valid'], 'scientifically_ready': False,
        'simulation_sha256': sha256(args.case / 'synthetic/simulation.json'),
        'estimate_sha256': sha256(args.output / 'estimate.json'), 'testing_data_access': False, 'final_data_access': False})


if __name__ == '__main__':
    main()
