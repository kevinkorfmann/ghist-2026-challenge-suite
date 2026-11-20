"""Report root's separate practical criteria without changing strict scientific gates."""
import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path

from ghist_lowcoverage.pipeline import sha256


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert args.root.name == 'lowcoverage-20261120' and not args.output.exists()
    policy = args.root / 'task14_pragmatic_usefulness_addendum_v1.json'
    assert sha256(policy) == '6318852e95ddc59004470a5e4971e3625fcb68bf60975b111a1191b272204f13'
    report_path = args.root / 'task14_pragmatic_development_gate_audit_v1.json'
    report = json.loads(report_path.read_text())
    assert report['freeze_sha256'] == '7a98a0e446cf5b5221fef48bace8743b5a604f7cd807fdff0c26df3df001b1bb'
    assert report['mode'] == 'development' and report['all_cohorts_completed']
    assert report['primary_target_count'] == 6
    assert set(report['cohorts']) == {'split_r0_p0_v1', 'split_r0_p1_v1'}
    keys = {'pop1_population_size', 'pop2_population_size', 'generations_since_split'}
    cases = {}
    for name, row in report['cohorts'].items():
        assert set(row['estimates']) == {'seed1', 'seed2', 'continued'}
        assert set(row['comparisons']) == {'seed2', 'continued'}
        valid = all(set(estimate) == keys and all(not isinstance(value, bool) and
                    math.isfinite(value) and value > 0 for value in estimate.values())
                    for estimate in row['estimates'].values())
        sensitivity = max(abs(value) for comparison in row['comparisons'].values()
                          for value in comparison['target_relative_change'].values())
        error = max(abs(value) for errors in row['relative_errors'].values() for value in errors.values())
        practical = {
            'valid_all_three_outputs': valid,
            'pre_post_provenance': row['gates']['pre_post_provenance'],
            'formatter_schema': row['gates']['formatter_schema'],
            'interior_success_and_two_optimizer_agreement': row['gates']['two_optimizer_agreement'],
            'exposure_at_least_99_percent': min(row['callable_fraction'].values()) >= .99,
            'all_target_errors_at_most_30_percent': error <= .30,
            'all_target_sensitivities_at_most_5_percent': sensitivity <= .05,
        }
        cases[name] = {'practical_criteria': practical, 'all_practical_criteria_passed': all(practical.values()),
            'maximum_absolute_target_error': error, 'maximum_absolute_target_sensitivity': sensitivity,
            'strict_failed_gates_retained': row['failed_gates'], 'strict_gates': row['gates'],
            'comparisons': row['comparisons']}
    result = {'created_at_utc': datetime.now(timezone.utc).isoformat(),
        'scope': 'Two previously inspected development cases; not fresh confirmation.',
        'freeze_sha256': report['freeze_sha256'], 'strict_report_sha256': sha256(report_path),
        'prospective_practical_addendum_sha256': sha256(policy), 'assessor_sha256': sha256(Path(__file__)),
        'cases': cases, 'all_practical_criteria_passed': all(row['all_practical_criteria_passed'] for row in cases.values()),
        'strict_development_gates_passed': report['all_registered_gates_passed'],
        'strict_qualification': False, 'scientifically_ready': False,
        'root_independent_audit_and_explicit_decision_required': True,
        'final_data_access': False, 'submitted': False}
    with args.output.open('x') as stream:
        stream.write(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
