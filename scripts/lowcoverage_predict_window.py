"""Explicit-input Task14 inference using the immutable v7 confirmation estimator.

No downloader, challenge truth reader, or submission transport. The coordinating
root agent owns any challenge-data execution after confirmation and audit.
"""
import argparse
import importlib.metadata as metadata
import json
from pathlib import Path
import sys

from ghist_lowcoverage.pipeline import audit_inputs, sha256, write_json
from ghist_lowcoverage.submission import format_submission_draft
from ghist_lowcoverage.winsfs import make_safs, run_window_comparison


FROZEN_SHA = '1f61090830f5d33cd62eb117488d32f906e670d8819add5f30edb1b26b0b11bc'


def verify_frozen(root):
    path = root / 'task14_confirmation_freeze_v1.json'
    if sha256(path) != FROZEN_SHA:
        raise ValueError('confirmation freeze checksum mismatch')
    frozen = json.loads(path.read_text())
    for name, expected in frozen['bindings'].items():
        if sha256(root / name) != expected:
            raise ValueError(f'frozen binding changed: {name}')
    for name, expected in frozen['source_files'].items():
        if sha256(root / 'source_v7' / name) != expected:
            raise ValueError(f'frozen source changed: {name}')
    runtime = json.loads((root / 'runtime_task14_confirmation_v1.json').read_text())
    for name, expected in runtime['files'].items():
        if sha256(Path(name)) != expected:
            raise ValueError(f'frozen runtime changed: {name}')
    if sys.version != runtime['python'] or any(metadata.version(name) != version
                                               for name, version in runtime['packages'].items()):
        raise ValueError('frozen runtime version mismatch')
    # Verify the modules actually imported, not just a dormant copy on disk.
    import ghist_lowcoverage.pipeline as pipeline
    import ghist_lowcoverage.submission as submission
    import ghist_lowcoverage.winsfs as winsfs
    for module in [pipeline, submission, winsfs]:
        relative = f'src/ghist_lowcoverage/{Path(module.__file__).name}'
        if Path(module.__file__).resolve() != (root / 'source_v7' / relative).resolve():
            raise ValueError('PYTHONPATH did not select the immutable source_v7 estimator')


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
    parser.add_argument('--directory', type=Path, required=True, help='BAMs, indices and popfile.txt')
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root, directory, reference = [path.resolve() for path in [args.frozen_root, args.directory, args.reference]]
    if args.output.exists():
        raise FileExistsError(args.output)
    verify_frozen(root)
    before = bind_inputs(directory, reference)
    args.output.mkdir(parents=True, exist_ok=False)
    write_json(args.output / 'initial_binding.json', {'freeze_sha256': FROZEN_SHA, 'input_binding': before,
        'wrapper_sha256': sha256(Path(__file__)), 'primary_estimate': 'seed1', 'threads': 4, 'max_epochs': 50})
    audit, indices = make_safs(directory, reference, 'split_isolation', args.output / 'safs',
        root / 'external/angsd/angsd', threads=4, expected_length=10_000_000)
    comparison = run_window_comparison(audit, indices, args.output / 'window_comparison',
        root / 'external/winsfs-release-v0.7.0/winsfs', threads=4, max_epochs=50)
    verify_frozen(root)
    after = bind_inputs(directory, reference)
    if before != after or audit != before['audit']:
        raise ValueError('input files changed during inference')
    gates = numerical_gates(comparison)
    result = {'freeze_sha256': FROZEN_SHA, 'primary_estimate': comparison['estimates']['seed1'],
        'numerical_gates': gates, 'all_numerical_gates_passed': all(gates.values()),
        'source_runtime_input_pre_post_verified': True, 'wrapper_sha256': sha256(Path(__file__)),
        'initial_binding_sha256': sha256(args.output / 'initial_binding.json'),
        'output_bindings': {name: sha256(args.output / name) for name in [
            'safs/inputs.json', 'safs/protocol.json', 'window_comparison/protocol.json', 'window_comparison/comparison.json']},
        'submitted': False, 'scientific_qualification': 'Requires separate all-six confirmation and independent audit.'}
    if all(gates.values()):
        draft = format_submission_draft('split_isolation', result['primary_estimate'],
            root / 'official_templates/split_isolation_low_coverage.submission.yaml')
        draft_path = args.output / 'task14.submission.draft.yaml'
        draft_path.write_text(draft)
        result['draft_sha256'] = sha256(draft_path)
    write_json(args.output / 'execution_binding.json', result)
    print(json.dumps({'all_numerical_gates_passed': all(gates.values()),
                      'failed_gates': [name for name, passed in gates.items() if not passed]}))


if __name__ == '__main__':
    main()
