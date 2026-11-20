"""Refit existing independent synthetic spectra with the corrected optimizer policy."""
import argparse
import json
from pathlib import Path

import numpy as np

from ghist_lowcoverage.pipeline import TASKS, fit_spectrum, sha256, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--sfs', type=Path, help='Existing inferred SFS; absent means synthetic genotype-truth control')
    parser.add_argument('--restarts', type=int, default=8)
    args = parser.parse_args()
    record_path = args.case / 'synthetic/simulation.json'
    simulation = json.loads(record_path.read_text())
    if simulation['testing_data_access'] or simulation['final_data_access']:
        raise ValueError('synthetic-only diagnostic')
    args.output.mkdir(parents=True, exist_ok=False)
    spectrum_path = args.sfs or args.case / 'synthetic/truth_folded_sfs.npy'
    if args.sfs:
        values = np.loadtxt(args.sfs)
    else:
        values = np.load(spectrum_path).astype(float)
        # tskit site AFS omits invariant sequence; restore measured simulation length.
        values.flat[0] = simulation['length'] - values.sum() + values.flat[0]
    estimate = fit_spectrum(values, TASKS[simulation['task']]['sample_sizes'], simulation['task'], restarts=args.restarts)
    write_json(args.output / 'estimate.json', estimate)
    errors = {k: estimate['targets'][k] / value - 1 for k, value in simulation['truth'].items()}
    write_json(args.output / 'scored.json', {'source_sha256': sha256(Path(__file__)),
        'spectrum_sha256': sha256(spectrum_path), 'simulation_sha256': sha256(record_path),
        'spectrum_kind': 'inferred' if args.sfs else 'synthetic_genotype_truth_control',
        'truth': simulation['truth'], 'estimate': estimate['targets'], 'relative_errors': errors,
        'optimizer_agreement': estimate['optimizer_agreement'], 'scientifically_ready': False,
        'testing_data_access': False, 'final_data_access': False})


if __name__ == '__main__':
    main()
