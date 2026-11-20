"""Evaluate the prospectively registered v7 full-width development gates.

This reads existing synthetic results only; it cannot promote a task to ready.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.root.name != 'lowcoverage-20261120' or args.output.exists():
        raise ValueError('isolated root and new output required')
    protocol = args.root / 'winsfs_full_development_v1.json'
    expected = '6bb7811d87d0028dfa5d817886e5d26dbe171585e7939b4fc027e8e28d4415f6'
    if digest(protocol) != expected:
        raise ValueError('these gates apply only to the registered development policy')
    report = {'protocol_sha256': expected, 'assessor_sha256': digest(Path(__file__)),
              'scientifically_ready': False, 'fresh_confirmation_completed': False, 'tasks': {}}
    for task in ['bottleneck', 'split_isolation']:
        directory = args.root / f'window_full_{task}_v7'
        scored = json.loads((directory / 'scored.json').read_text())
        binding = json.loads((directory / 'execution_binding.json').read_text())
        comparison = scored['comparison']
        assert binding['scored_result_sha256'] == digest(directory / 'scored.json')
        assert binding['protocol_sha256'] == expected
        errors = {label: {key: estimate[key] / truth - 1 for key, truth in scored['truth'].items()}
                  for label, estimate in comparison['estimates'].items()}
        assert errors == scored['relative_errors']
        gates = {
            'pre_post_binding': binding['source_and_input_pre_post_verification_passed'],
            'exposure': all(value >= .99 for value in comparison['callable_fraction'].values()),
            'observed_stopping': all(comparison['early_stop_observed'].values()),
            'two_optimizer_agreement': all(value['passed'] for value in comparison['optimizer_agreement'].values()),
            'target_stability': all(abs(value) < .01 for row in comparison['comparisons'].values()
                                    for value in row['target_relative_change'].values()),
            'frequency_L1': all(row['frequency_L1'] < .001 for row in comparison['comparisons'].values()),
            'segregating_relative_L1': all(row['segregating_relative_L1'] < .05 for row in comparison['comparisons'].values()),
            'parameter_recovery': all(abs(value) <= .30 for row in errors.values() for value in row.values()),
        }
        report['tasks'][task] = {
            'gates': gates, 'full_bam_development_gates_passed': all(gates.values()),
            'failed_gates': [name for name, passed in gates.items() if not passed],
            'estimates': comparison['estimates'], 'relative_errors': errors,
            'comparisons': comparison['comparisons'], 'callable_fraction': comparison['callable_fraction'],
            'documented_demography_root_sum_squared_relative_errors': {
                label: math.sqrt(sum(value * value for value in row.values())) for label, row in errors.items()},
            'metric_applicability': 'Published demographic-history formula; application to added Tasks13/14 inferred, not separately confirmed.',
            'scored_sha256': digest(directory / 'scored.json'),
            'execution_binding_sha256': digest(directory / 'execution_binding.json'),
            'scientifically_ready': False,
        }
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps({task: row['failed_gates'] for task, row in report['tasks'].items()}))


if __name__ == '__main__':
    main()
