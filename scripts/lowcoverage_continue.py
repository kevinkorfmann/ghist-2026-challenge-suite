"""Continue original smoke realSFS fit and independently cross-check optimization."""
import argparse
import json
from pathlib import Path
import subprocess

import numpy as np

from ghist_lowcoverage.pipeline import fit_spectrum, realsfs_diagnostics, sha256, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--iterations', type=int, default=5000)
    parser.add_argument('--restarts', type=int, default=8)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    original = args.case / 'inference'
    audit = json.loads((original / 'inputs.json').read_text())
    for item in list(audit['samples'].values()) + [{'path': audit['reference'], 'sha256': audit['reference_sha256']}]:
        if sha256(Path(item['path'])) != item['sha256']:
            raise ValueError('original inference input changed')
    command = json.loads((original / 'protocol.json').read_text())['commands'][-1]
    command[command.index('-maxIter') + 1] = str(args.iterations)
    command += ['-sfs', str(original / 'folded.sfs')]
    write_json(args.output / 'protocol.json', {'command': command, 'source_sha256': sha256(Path(__file__)),
        'original_sfs_sha256': sha256(original / 'folded.sfs'), 'original_input_audit_sha256': sha256(original / 'inputs.json'),
        'case': str(args.case), 'testing_data_access': False, 'final_data_access': False})
    with (args.output / 'folded.sfs').open('w') as out, (args.output / 'realsfs.log').open('w') as log:
        subprocess.run(command, stdout=out, stderr=log, check=True)
    before = np.loadtxt(original / 'folded.sfs')
    after = np.loadtxt(args.output / 'folded.sfs')
    estimates = {}
    for name, spectrum in [('initial', before), ('continued', after)]:
        estimates[name] = fit_spectrum(spectrum, audit['sample_sizes'], audit['task'], restarts=args.restarts)
        write_json(args.output / f'{name}_estimate.json', estimates[name])
    truth = json.loads((args.case / 'synthetic/simulation.json').read_text())['truth']
    errors = {name: {key: result['targets'][key] / true - 1 for key, true in truth.items()}
              for name, result in estimates.items()}
    relative_change = {key: estimates['continued']['targets'][key] / estimates['initial']['targets'][key] - 1 for key in truth}
    diagnostics = realsfs_diagnostics((args.output / 'realsfs.log').read_text())
    write_json(args.output / 'result.json', {'task': audit['task'], 'truth': truth, 'relative_errors': errors,
        'target_change_after_continuation': relative_change, 'spectrum_count_L1_change': float(np.abs(after-before).sum()),
        'realsfs_diagnostics': diagnostics, 'optimizer_agreement': estimates['continued']['optimizer_agreement'],
        'numerically_valid': diagnostics['convergence_observed'] and estimates['continued']['optimizer_agreement']['passed'],
        'scientifically_ready': False, 'official_score': False})


if __name__ == '__main__':
    main()
