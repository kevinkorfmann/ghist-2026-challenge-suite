"""Recompute fixed stabilization gates for both development or all six fresh cases."""
import argparse
import json
from pathlib import Path
from statistics import median

from ghist_lowcoverage.pipeline import sha256, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--freeze-sha', required=True)
    parser.add_argument('--mode', choices=['development', 'confirmation'], required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert args.root.name == 'lowcoverage-20261120' and not args.output.exists()
    assert sha256(args.root / 'task14_extension_freeze_v2.json') == args.freeze_sha
    plan_path = args.root / 'task14_extension_protocol_v2.json'
    assert sha256(plan_path) == '861c6b5618c2afd665c72d8dd4773c16e592cac587b410267f9bd32cf287bef2'
    plan = json.loads(plan_path.read_text())
    regimes = range(1) if args.mode == 'development' else range(3)
    report = {'freeze_sha256': args.freeze_sha, 'assessor_sha256': sha256(Path(__file__)),
        'protocol_sha256': sha256(plan_path), 'mode': args.mode, 'cohorts': {}, 'scientifically_ready': False,
        'independent_audit_pending': True, 'previous_confirmation_status': 'Failed and retained unchanged.'}
    primary, strata = [], {0: [], 1: []}
    for regime in regimes:
        for replicate in range(2):
            stem = f'split_r{regime}_p{replicate}_v1'
            prefix = 'full_confirmation_' if args.mode == 'development' else 'full_extended_confirmation_'
            case = args.root / f'{prefix}{stem}'
            directory = args.root / f'window_extended_{args.mode}_{stem}'
            files = [directory / name for name in ['initial_binding.json', 'execution_binding.json', 'scored.json']]
            if any(not path.is_file() for path in files):
                report['cohorts'][stem] = {'all_gates_passed': False,
                    'missing_artifacts': [str(path.relative_to(args.root)) for path in files if not path.is_file()]}
                continue
            initial, executed, scored = [json.loads(path.read_text()) for path in files]
            assert initial['freeze_sha256'] == executed['freeze_sha256'] == scored['freeze_sha256'] == args.freeze_sha
            assert executed['initial_binding_sha256'] == sha256(files[0])
            assert initial['input_binding'] == executed['input_binding']
            for name, expected in executed['output_bindings'].items():
                assert sha256(directory / name) == expected
            simulation_path = case / 'synthetic/simulation.json'
            assert sha256(simulation_path) == scored['simulation_sha256'] == initial['input_binding']['simulation_sha256']
            simulation = json.loads(simulation_path.read_text())
            assert not simulation['testing_data_access'] and not simulation['final_data_access']
            assert simulation['length'] == 10_000_000 and simulation['truth'] == scored['truth']
            assert initial['input_binding']['audit'] == json.loads((directory / 'safs/inputs.json').read_text())
            if args.mode == 'confirmation':
                prepared = json.loads((case / 'prepare_binding.json').read_text())
                assert prepared['freeze_sha256'] == args.freeze_sha and prepared['input_binding'] == initial['input_binding']
                expected_design = dict(plan['confirmation_designs']['split_isolation'][regime])
                expected_design['seed'] = expected_design.pop('seeds')[replicate]
                assert simulation['design'] == expected_design
            c = scored['comparison']
            labels = {'seed1', 'seed2', 'continued'}
            assert all(set(c[name]) == labels for name in ['estimates', 'callable_fraction', 'early_stop_observed', 'optimizer_agreement'])
            assert set(c['comparisons']) == {'seed2', 'continued'}
            errors = {label: {key: estimate[key] / truth - 1 for key, truth in scored['truth'].items()}
                      for label, estimate in c['estimates'].items()}
            assert errors == scored['relative_errors']
            gates = {
                'pre_post_provenance': executed['source_runtime_input_pre_post_verified'],
                'formatter_schema': executed['synthetic_formatter_schema_passed'],
                'exposure': all(value >= .99 for value in c['callable_fraction'].values()),
                'observed_stopping': all(c['early_stop_observed'].values()),
                'two_optimizer_agreement': all(row['passed'] for row in c['optimizer_agreement'].values()),
                'target_stability': all(abs(value) < .01 for row in c['comparisons'].values()
                                        for value in row['target_relative_change'].values()),
                'frequency_L1': all(row['frequency_L1'] < .001 for row in c['comparisons'].values()),
                'segregating_relative_L1': all(row['segregating_relative_L1'] < .05 for row in c['comparisons'].values()),
                'primary_parameter_recovery': all(abs(value) <= .30 for value in errors['seed1'].values()),
                'all_estimator_runs_parameter_recovery': all(abs(value) <= .30 for row in errors.values() for value in row.values()),
            }
            values = [abs(value) for value in errors['seed1'].values()]
            assert len(values) == 3
            primary.extend(values)
            strata[replicate].extend(values)
            report['cohorts'][stem] = {'gates': gates, 'all_gates_passed': all(gates.values()),
                'failed_gates': [name for name, passed in gates.items() if not passed], 'design': simulation['design'],
                'noise_profile': {'quality': simulation['quality'], 'depth_sigma': simulation['individual_depth_lognormal_sigma']},
                'estimates': c['estimates'], 'relative_errors': errors, 'comparisons': c['comparisons'],
                'callable_fraction': c['callable_fraction'],
                'evidence_sha256': {str(path.relative_to(args.root)): sha256(path) for path in files}}
    expected_count = 6 if args.mode == 'development' else 18
    report['primary_target_count'] = len(primary)
    report['all_cohorts_completed'] = len(primary) == expected_count
    report['primary_median_absolute_relative_error'] = median(primary) if len(primary) == expected_count else None
    report['stratum_medians'] = {f'replicate{key}': median(values) if len(values) == expected_count // 2 else None
                                 for key, values in strata.items()}
    report['median_gate_required'] = args.mode == 'confirmation'
    report['median_recovery_gate_passed'] = len(primary) == expected_count and (
        args.mode == 'development' or median(primary) <= .20)
    report['all_registered_gates_passed'] = report['median_recovery_gate_passed'] and all(
        row['all_gates_passed'] for row in report['cohorts'].values())
    write_json(args.output, report)
    print(json.dumps({'all_registered_gates_passed': report['all_registered_gates_passed'],
        'primary_target_count': len(primary), 'median': report['primary_median_absolute_relative_error'],
        'failures': {name: row.get('failed_gates', row.get('missing_artifacts'))
                     for name, row in report['cohorts'].items() if not row['all_gates_passed']}}))


if __name__ == '__main__':
    main()
