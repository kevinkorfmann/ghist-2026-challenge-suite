"""Assess all six frozen Task14 synthetic cohorts without promoting readiness."""
import argparse
import hashlib
import json
from pathlib import Path
from statistics import median


FROZEN_SHA = '1f61090830f5d33cd62eb117488d32f906e670d8819add5f30edb1b26b0b11bc'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.root.name != 'lowcoverage-20261120' or args.output.exists():
        raise ValueError('isolated experiment root and new output required')
    frozen_path = args.root / 'task14_confirmation_freeze_v1.json'
    assert sha(frozen_path) == FROZEN_SHA
    frozen = json.loads(frozen_path.read_text())
    report = {'freeze_sha256': FROZEN_SHA, 'assessor_sha256': sha(Path(__file__)), 'cohorts': {},
              'scientifically_ready': False, 'independent_final_audit_pending': True,
              'scope': 'All six prespecified synthetic Task14 cohorts; no challenge data.'}
    primary, strata = [], {0: [], 1: []}
    for regime in range(3):
        for replicate in range(2):
            stem = f'split_r{regime}_p{replicate}_v1'
            case = args.root / f'full_confirmation_{stem}'
            output = args.root / f'window_confirmation_{stem}'
            required = [case / 'synthetic/simulation.json', case / 'prepare_binding.json',
                        args.root / f'full_confirmation_{stem}.initial.json',
                        args.root / f'window_confirmation_{stem}.initial.json',
                        output / 'execution_binding.json', output / 'scored.json']
            missing = [str(path.relative_to(args.root)) for path in required if not path.is_file()]
            if missing:
                report['cohorts'][stem] = {'all_gates_passed': False, 'missing_artifacts': missing}
                continue
            simulation, prepared, prep_initial, infer_initial, executed, scored = [
                json.loads(path.read_text()) for path in required]
            design = dict(frozen['confirmation_designs'][regime])
            design['seed'] = design.pop('seeds')[replicate]
            assert simulation['design'] == design and simulation['length'] == 10_000_000
            noise = frozen['noise_profiles'][f'replicate{replicate}']
            assert simulation['quality'] == noise['quality']
            assert simulation['individual_depth_lognormal_sigma'] == noise['depth_sigma']
            assert not simulation['testing_data_access'] and not simulation['final_data_access']
            for record in [prepared, prep_initial, infer_initial, executed]:
                assert record['freeze_sha256'] == FROZEN_SHA
                assert record['runtime_source_tools_formatter_verified']
                assert record['regime'] == regime and record['replicate'] == replicate
            assert prepared['input_binding'] == infer_initial['input_binding'] == executed['input_binding']
            assert prepared['input_binding']['simulation_sha256'] == sha(required[0])
            assert prepared['input_binding']['audit'] == json.loads((output / 'safs/inputs.json').read_text())
            for name, expected in executed['output_bindings'].items():
                assert sha(output / name) == expected
            assert scored['truth'] == simulation['truth']
            comparison = scored['comparison']
            assert set(comparison['estimates']) == {'seed1', 'seed2', 'continued'}
            errors = {label: {key: estimate[key] / truth - 1 for key, truth in scored['truth'].items()}
                      for label, estimate in comparison['estimates'].items()}
            assert errors == scored['relative_errors']
            gates = {
                'provenance_bindings': True,
                'formatter_schema': executed['synthetic_formatter_schema_passed'],
                'exposure': all(value >= .99 for value in comparison['callable_fraction'].values()),
                'observed_stopping': all(comparison['early_stop_observed'].values()),
                'two_optimizer_agreement': all(value['passed'] for value in comparison['optimizer_agreement'].values()),
                'target_stability': all(abs(value) < .01 for row in comparison['comparisons'].values()
                                        for value in row['target_relative_change'].values()),
                'frequency_L1': all(row['frequency_L1'] < .001 for row in comparison['comparisons'].values()),
                'segregating_relative_L1': all(row['segregating_relative_L1'] < .05 for row in comparison['comparisons'].values()),
                'primary_parameter_recovery': all(abs(value) <= .30 for value in errors['seed1'].values()),
                'all_estimator_runs_parameter_recovery': all(abs(value) <= .30 for row in errors.values() for value in row.values()),
            }
            absolute = [abs(value) for value in errors['seed1'].values()]
            assert len(absolute) == 3
            primary.extend(absolute)
            strata[replicate].extend(absolute)
            report['cohorts'][stem] = {'gates': gates, 'all_gates_passed': all(gates.values()),
                'failed_gates': [key for key, passed in gates.items() if not passed],
                'design': design, 'noise_profile': noise, 'estimates': comparison['estimates'],
                'relative_errors': errors, 'comparisons': comparison['comparisons'],
                'callable_fraction': comparison['callable_fraction'],
                'evidence_sha256': {str(path.relative_to(args.root)): sha(path) for path in required}}
    report['primary_target_count'] = len(primary)
    report['primary_median_absolute_relative_error'] = median(primary) if len(primary) == 18 else None
    report['stratum_medians'] = {f'replicate{key}': median(values) if len(values) == 9 else None
                                 for key, values in strata.items()}
    report['all_six_cohorts_completed'] = len(primary) == 18
    report['median_recovery_gate_passed'] = len(primary) == 18 and median(primary) <= .20
    report['all_registered_gates_passed'] = report['median_recovery_gate_passed'] and all(
        row['all_gates_passed'] for row in report['cohorts'].values())
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'all_registered_gates_passed': report['all_registered_gates_passed'],
                      'primary_target_count': len(primary),
                      'primary_median_absolute_relative_error': report['primary_median_absolute_relative_error'],
                      'failures': {name: row.get('failed_gates', row.get('missing_artifacts'))
                                   for name, row in report['cohorts'].items() if not row['all_gates_passed']}}))


if __name__ == '__main__':
    main()
