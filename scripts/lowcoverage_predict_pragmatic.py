"""Root-only explicit-input pragmatic Task14 wrapper; all strict failures remain reported.

The frozen engine is unchanged. A separate hash-bound root decision is required
before input access. This module never downloads, reads truth, or submits.
"""
import argparse
import json
from pathlib import Path

from ghist_lowcoverage.pipeline import audit_inputs, sha256, write_json
from ghist_lowcoverage.submission import format_submission_draft
from ghist_lowcoverage.winsfs import make_safs
from ghist_lowcoverage.stabilized import run_stabilized_comparison
from scripts.lowcoverage_pragmatic_guard import check_clock, verify


def bind_inputs(directory, reference):
    audit = audit_inputs(directory, reference, 'split_isolation', 10_000_000)
    indices = [Path(str(reference) + '.fai')]
    for row in audit['samples'].values():
        bam = Path(row['path'])
        candidates = [Path(str(bam) + '.bai'), bam.with_suffix('.bai'), Path(str(bam) + '.csi')]
        existing = [path for path in candidates if path.is_file()]
        if not existing:
            raise ValueError(f'cannot bind BAM index: {bam.name}')
        indices.extend(existing)
    return {'audit': audit, 'index_sha256': {str(path): sha256(path) for path in indices}}


def numerical_gates(comparison):
    labels = {'seed1', 'seed2', 'continued'}
    if any(set(comparison[name]) != labels for name in [
        'estimates', 'callable_fraction', 'early_stop_observed', 'optimizer_agreement']):
        raise ValueError('all three frozen estimator runs are required')
    if set(comparison['comparisons']) != {'seed2', 'continued'}:
        raise ValueError('both frozen estimator comparisons are required')
    targets = {'pop1_population_size', 'pop2_population_size', 'generations_since_split'}
    if any(set(row['target_relative_change']) != targets for row in comparison['comparisons'].values()):
        raise ValueError('all three physical targets must be compared')
    return {
        'exposure': all(value >= .99 for value in comparison['callable_fraction'].values()),
        'observed_stopping': all(comparison['early_stop_observed'].values()),
        'two_optimizer_agreement': all(row['passed'] for row in comparison['optimizer_agreement'].values()),
        'target_stability': all(abs(value) < .01 for row in comparison['comparisons'].values()
                               for value in row['target_relative_change'].values()),
        'frequency_L1': all(row['frequency_L1'] < .001 for row in comparison['comparisons'].values()),
        'segregating_relative_L1': all(row['segregating_relative_L1'] < .05 for row in comparison['comparisons'].values()),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--frozen-root', type=Path, required=True)
    parser.add_argument('--freeze-sha', required=True)
    parser.add_argument('--root-decision', type=Path, required=True)
    parser.add_argument('--decision-sha', required=True)
    parser.add_argument('--directory', type=Path, required=True)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root, directory, reference = [path.resolve() for path in [args.frozen_root, args.directory, args.reference]]
    if args.output.exists():
        raise FileExistsError(args.output)
    clock_record = check_clock(root, 'final')
    verify(root, args.freeze_sha)
    if sha256(args.root_decision) != args.decision_sha:
        raise ValueError('root decision hash mismatch')
    decision = json.loads(args.root_decision.read_text())
    report_path = root / 'task14_pragmatic_development_gate_audit_v1.json'
    report = json.loads(report_path.read_text())
    if not (decision['authorize_final_inference'] is True and decision['strict_qualification'] is False
            and decision['freeze_sha256'] == args.freeze_sha
            and decision['wrapper_sha256'] == sha256(Path(__file__))
            and decision['development_report_sha256'] == sha256(report_path)
            and decision['primary_estimate'] == 'seed1'
            and report['freeze_sha256'] == args.freeze_sha
            and report['mode'] == 'development' and report['all_cohorts_completed']
            and report['primary_target_count'] == 6):
        raise ValueError('explicit root pragmatic decision after both complete development cases is required')
    before = bind_inputs(directory, reference)
    args.output.mkdir(parents=True, exist_ok=False)
    write_json(args.output / 'initial_binding.json', {'freeze_sha256': args.freeze_sha,
        'input_binding': before, 'wrapper_sha256': sha256(Path(__file__)),
        'primary_estimate': 'seed1', 'threads': 4, 'max_epochs': 500,
        'development_report_sha256': sha256(report_path), 'root_decision_sha256': args.decision_sha,
        'clock_guard': clock_record, 'strict_qualification': False})
    audit, indices = make_safs(directory, reference, 'split_isolation', args.output / 'safs',
        root / 'external/angsd/angsd', threads=4, expected_length=10_000_000)
    comparison = run_stabilized_comparison(audit, indices, args.output / 'window_comparison',
        root / 'external/winsfs-release-v0.7.0/winsfs', threads=4, max_epochs=500)
    verify(root, args.freeze_sha)
    initial = json.loads((args.output / 'initial_binding.json').read_text())
    if sha256(report_path) != initial['development_report_sha256'] or sha256(args.root_decision) != args.decision_sha:
        raise ValueError('development evidence or root decision changed during inference')
    after = bind_inputs(directory, reference)
    if before != after or audit != before['audit']:
        raise ValueError('input files changed during inference')
    for row in json.loads((args.output / 'window_comparison/protocol.json').read_text())['saf_files']:
        if sha256(Path(row['path'])) != row['sha256']:
            raise ValueError('SAF input changed during inference')
    gates = numerical_gates(comparison)
    result = {'freeze_sha256': args.freeze_sha, 'primary_estimate': comparison['estimates']['seed1'],
        'numerical_gates': gates, 'all_numerical_gates_passed': all(gates.values()),
        'source_runtime_input_pre_post_verified': True, 'wrapper_sha256': sha256(Path(__file__)),
        'initial_binding_sha256': sha256(args.output / 'initial_binding.json'),
        'root_decision_sha256': args.decision_sha,
        'output_bindings': {name: sha256(args.output / name) for name in [
            'safs/inputs.json', 'safs/protocol.json', 'window_comparison/protocol.json', 'window_comparison/comparison.json']},
        'submitted': False, 'strict_qualification': False, 'scientifically_ready': False,
        'scientific_qualification': 'Pragmatic candidate only; original six confirmation and cap200 failures preserved; no fresh cap500 confirmation.',
        'submission_requires_root_review_of_actual_numerical_gates': True}
    draft = format_submission_draft('split_isolation', result['primary_estimate'],
        root / 'official_templates/split_isolation_low_coverage.submission.yaml')
    draft_path = args.output / 'task14.pragmatic.submission.draft.yaml'
    with draft_path.open('x') as stream:
        stream.write(draft)
    result['draft_sha256'] = sha256(draft_path)
    write_json(args.output / 'execution_binding.json', result)
    print(json.dumps({'all_numerical_gates_passed': all(gates.values()),
        'failed_gates': [name for name, passed in gates.items() if not passed],
        'strict_qualification': False, 'submitted': False}))


if __name__ == '__main__':
    main()
