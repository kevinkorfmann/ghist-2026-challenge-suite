"""Assess the preregistered Task13 full-batch development diagnostic."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


REGISTERED_SHA = 'd25f94eae94fd168c6bb845f295ceb6eba2b97d7af10e16ad669c0f9667c7120'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compare(a, b, left, right):
    distances = {
        'frequency_L1': float(np.abs(a / a.sum() - b / b.sum()).sum()),
        'segregating_relative_L1': float(np.abs(a[1:] - b[1:]).sum() / a[1:].sum()),
        'target_relative_change': {key: right['targets'][key] / value - 1
                                   for key, value in left['targets'].items()}}
    distances['all_stability_gates_passed'] = (
        distances['frequency_L1'] < 1e-5 and distances['segregating_relative_L1'] < .01
        and all(abs(value) < .01 for value in distances['target_relative_change'].values()))
    return distances


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.root.name != 'lowcoverage-20261120' or args.output.exists():
        raise ValueError('isolated experiment root and new output required')
    registered = args.root / 'realsfs_refinement_development_v1.json'
    assert sha(registered) == REGISTERED_SHA
    report = {'registered_protocol_sha256': REGISTERED_SHA, 'assessor_sha256': sha(Path(__file__)),
              'scope': 'Previously inspected development only; prior v7/truth-control failures remain unchanged.',
              'scientifically_ready': False, 'initializations': {}}
    estimates, spectra = {}, {}
    for label in ['seed1', 'seed2', 'continued']:
        directory = args.root / f'full_refinement_bottleneck_{label}_v1'
        required = [directory / 'scored.json', directory / 'refined.estimate.json',
                    directory / 'continued.estimate.json', directory / 'refined.sfs', directory / 'continued.sfs']
        missing = [str(path.relative_to(args.root)) for path in required if not path.is_file()]
        if missing:
            report['initializations'][label] = {'all_gates_passed': False, 'missing_artifacts': missing}
            continue
        scored = json.loads(required[0].read_text())
        phase_estimates = {phase: json.loads((directory / f'{phase}.estimate.json').read_text())
                           for phase in ['refined', 'continued']}
        phase_spectra = {phase: np.loadtxt(directory / f'{phase}.sfs') for phase in phase_estimates}
        assert scored['protocol_sha256'] == sha(directory / 'protocol.json')
        assert scored['initial_binding_sha256'] == sha(directory / 'initial_binding.json')
        assert json.loads((directory / 'protocol.json').read_text())['registered_protocol_sha256'] == REGISTERED_SHA
        errors = {phase: {key: estimate['targets'][key] / truth - 1 for key, truth in scored['truth'].items()}
                  for phase, estimate in phase_estimates.items()}
        assert errors == scored['relative_errors']
        comparison = compare(phase_spectra['refined'], phase_spectra['continued'],
                             phase_estimates['refined'], phase_estimates['continued'])
        gates = {
            'pre_post_provenance': scored['source_and_input_pre_post_verification_passed'],
            'exposure': all(estimate['callable_fraction'] >= .99 for estimate in phase_estimates.values()),
            'observed_stopping_before_cap': all(estimate['realsfs_diagnostics']['convergence_observed']
                and estimate['realsfs_diagnostics']['last_iteration'] is not None
                and estimate['realsfs_diagnostics']['last_iteration'] < 250 for estimate in phase_estimates.values()),
            'two_optimizer_agreement': all(estimate['optimizer_agreement']['passed'] for estimate in phase_estimates.values()),
            'parameter_recovery': all(abs(value) <= .30 for row in errors.values() for value in row.values()),
            'continuation_stability': comparison['all_stability_gates_passed'],
        }
        report['initializations'][label] = {'gates': gates, 'all_gates_passed': all(gates.values()),
            'failed_gates': [name for name, passed in gates.items() if not passed],
            'estimates': {phase: estimate['targets'] for phase, estimate in phase_estimates.items()},
            'relative_errors': errors, 'continuation_comparison': comparison,
            'stopping_diagnostics': {phase: estimate['realsfs_diagnostics'] for phase, estimate in phase_estimates.items()},
            'evidence_sha256': {str(path.relative_to(args.root)): sha(path) for path in required}}
        estimates[label], spectra[label] = phase_estimates, phase_spectra
    report['all_initializations_completed'] = len(estimates) == 3
    report['between_initializations'] = {}
    if len(estimates) == 3:
        # The prospective protocol explicitly uses refined seed1 as the common anchor.
        report['between_initializations'] = {label: compare(spectra['seed1']['refined'], spectra[label]['continued'],
            estimates['seed1']['refined'], estimates[label]['continued']) for label in estimates}
    report['all_registered_diagnostic_gates_passed'] = len(estimates) == 3 and all(
        row['all_gates_passed'] for row in report['initializations'].values()) and all(
        row['all_stability_gates_passed'] for row in report['between_initializations'].values())
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'all_registered_diagnostic_gates_passed': report['all_registered_diagnostic_gates_passed'],
        'failures': {name: row.get('failed_gates', row.get('missing_artifacts'))
                     for name, row in report['initializations'].items() if not row['all_gates_passed']}}))


if __name__ == '__main__':
    main()
