"""Matched development comparison using existing synthetic GL2 SAF files."""
import argparse
import json
from pathlib import Path
import subprocess

import numpy as np

from ghist_lowcoverage.pipeline import fit_spectrum, realsfs_diagnostics, sha256, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', type=Path, required=True)
    parser.add_argument('--saf-directory', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--realsfs', type=Path, required=True)
    parser.add_argument('--threads', type=int, default=4)
    args = parser.parse_args()
    simulation = json.loads((args.case/'synthetic/simulation.json').read_text())
    if simulation['testing_data_access'] or simulation['final_data_access']:
        raise ValueError('independent synthetic comparison only')
    audit = json.loads((args.saf_directory/'inputs.json').read_text())
    indices = [args.saf_directory/f'population_{i}.saf.idx' for i in range(len(audit['population_order']))]
    args.output.mkdir(parents=True, exist_ok=False)
    command = [str(args.realsfs), *map(str, indices), '-fold', '1', '-P', str(args.threads),
               '-maxIter', '1000', '-tole', '1e-8', '-seed', '20261120']
    write_json(args.output/'protocol.json', {'source_sha256': sha256(Path(__file__)), 'command': command,
        'input_protocol_sha256': sha256(args.saf_directory/'protocol.json'),
        'realsfs_sha256': sha256(args.realsfs),
        'saf_files': [{'path': str(member), 'sha256': sha256(member)} for index in indices
                      for member in [index, Path(str(index).removesuffix('.idx')+'.gz'),
                                     Path(str(index).removesuffix('.idx')+'.pos.gz')]],
        'testing_data_access': False, 'final_data_access': False})
    with (args.output/'folded.sfs').open('w') as out, (args.output/'realsfs.log').open('w') as log:
        subprocess.run(command, stdout=out, stderr=log, check=True)
    estimate = fit_spectrum(np.loadtxt(args.output/'folded.sfs'), audit['sample_sizes'], audit['task'])
    estimate['realsfs_diagnostics'] = realsfs_diagnostics((args.output/'realsfs.log').read_text())
    write_json(args.output/'estimate.json', estimate)
    write_json(args.output/'scored.json', {'truth': simulation['truth'], 'estimate': estimate['targets'],
        'relative_errors': {k: estimate['targets'][k]/v-1 for k,v in simulation['truth'].items()},
        'optimizer_agreement': estimate['optimizer_agreement'], 'realsfs_diagnostics': estimate['realsfs_diagnostics'],
        'scientifically_ready': False, 'testing_data_access': False, 'final_data_access': False})


if __name__ == '__main__':
    main()
