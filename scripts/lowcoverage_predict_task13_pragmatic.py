"""Root-only Task13 fallback using a fixed, previously evaluated continuation.

No downloads, truth access, synthetic fitting, or submission. Development-based
method selection is explicit; neither recovery nor scientific readiness passed.
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

from ghist_lowcoverage.pipeline import audit_inputs, sha256, write_json
from ghist_lowcoverage.submission import format_submission_draft
from ghist_lowcoverage.winsfs import make_safs, run_window_comparison
from scripts.lowcoverage_pragmatic_guard import verify


FREEZE = '7a98a0e446cf5b5221fef48bace8743b5a604f7cd807fdff0c26df3df001b1bb'
TEMPLATE_SHA = 'd2d4ad595d4de0466fc4291d3f6550054da98676e42c662b6b57aaeabf5a9936'
PRIMARY = 'continued'


def bind_inputs(directory, reference):
    audit = audit_inputs(directory, reference, 'bottleneck', 10_000_000)
    indices = [Path(str(reference) + '.fai')]
    for row in audit['samples'].values():
        bam = Path(row['path'])
        candidates = [Path(str(bam) + '.bai'), bam.with_suffix('.bai'), Path(str(bam) + '.csi')]
        existing = [path for path in candidates if path.is_file()]
        if not existing:
            raise ValueError(f'cannot bind BAM index: {bam.name}')
        indices.extend(existing)
    return {'audit': audit, 'index_sha256': {str(path): sha256(path) for path in indices}}


def check_decision(decision, wrapper_sha):
    if not (decision['authorize_final_inference'] is True
            and decision['task'] == 'bottleneck'
            and decision['strict_qualification'] is False
            and decision['freeze_sha256'] == FREEZE
            and decision['wrapper_sha256'] == wrapper_sha
            and decision['primary_estimate'] == PRIMARY
            and decision['development_selection_disclosed'] is True):
        raise ValueError('explicit Task13 fallback decision is required')


def numerical_gates(comparison):
    labels = {'seed1', 'seed2', 'continued'}
    if any(set(comparison[name]) != labels for name in [
            'estimates', 'callable_fraction', 'early_stop_observed', 'optimizer_agreement']):
        raise ValueError('all three fixed estimator runs required')
    return {
        'exposure': all(value >= .99 for value in comparison['callable_fraction'].values()),
        'observed_stopping': all(comparison['early_stop_observed'].values()),
        'two_optimizer_agreement': all(value['passed'] for value in comparison['optimizer_agreement'].values()),
        'target_stability': all(abs(value) < .01 for row in comparison['comparisons'].values()
                                for value in row['target_relative_change'].values()),
        'frequency_L1': all(row['frequency_L1'] < .001 for row in comparison['comparisons'].values()),
        'segregating_relative_L1': all(row['segregating_relative_L1'] < .05 for row in comparison['comparisons'].values()),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['frozen-root', 'root-decision', 'directory', 'reference', 'output']:
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--decision-sha', required=True)
    args = parser.parse_args()
    root, directory, reference = [p.resolve() for p in [args.frozen_root, args.directory, args.reference]]
    if args.output.exists():
        raise FileExistsError(args.output)
    now = datetime.now(timezone.utc)
    if now > datetime.fromisoformat('2026-11-20T22:45:00+00:00'):
        raise RuntimeError('55 minute fallback must start by 22:45 UTC')
    verify(root, FREEZE)
    template = root / 'official_templates/bottleneck.submission.yaml'
    if sha256(template) != TEMPLATE_SHA or sha256(args.root_decision) != args.decision_sha:
        raise ValueError('template or root decision hash mismatch')
    wrapper_sha = sha256(Path(__file__))
    check_decision(json.loads(args.root_decision.read_text()), wrapper_sha)
    before = bind_inputs(directory, reference)
    args.output.mkdir(parents=True, exist_ok=False)
    write_json(args.output / 'initial_binding.json', {
        'freeze_sha256': FREEZE, 'input_binding': before, 'wrapper_sha256': wrapper_sha,
        'primary_estimate': PRIMARY, 'root_decision_sha256': args.decision_sha,
        'actual_start_utc': now.isoformat(), 'strict_qualification': False,
        'development_selection': 'Continued winsfs has RSS error .414327 versus .517388 for continued realSFS; no bias correction.',
        'known_development_failures': 'Continued decline fraction +41.21%; initial seed disagreement; no fresh confirmation.'})
    audit, indices = make_safs(directory, reference, 'bottleneck', args.output / 'safs',
        root / 'external/angsd/angsd', threads=4, expected_length=10_000_000)
    comparison = run_window_comparison(audit, indices, args.output / 'window_comparison',
        root / 'external/winsfs-release-v0.7.0/winsfs', threads=4, max_epochs=50)
    verify(root, FREEZE)
    if (sha256(template) != TEMPLATE_SHA or sha256(args.root_decision) != args.decision_sha
            or sha256(Path(__file__)) != wrapper_sha):
        raise ValueError('template or root decision changed during inference')
    if before != bind_inputs(directory, reference) or audit != before['audit']:
        raise ValueError('input files changed during inference')
    protocol = json.loads((args.output / 'window_comparison/protocol.json').read_text())
    for row in protocol['saf_files']:
        if sha256(Path(row['path'])) != row['sha256']:
            raise ValueError('SAF changed during inference')
    primary = comparison['estimates'][PRIMARY]
    draft = format_submission_draft('bottleneck', primary, template)
    draft_path = args.output / 'task13.pragmatic.submission.draft.yaml'
    with draft_path.open('x') as stream:
        stream.write(draft)
    gates = numerical_gates(comparison)
    result = {'task': 'bottleneck', 'freeze_sha256': FREEZE, 'primary_label': PRIMARY,
        'primary_estimate': primary, 'numerical_gates': gates,
        'all_numerical_gates_passed': all(gates.values()),
        'source_runtime_input_pre_post_verified': True, 'wrapper_sha256': wrapper_sha,
        'initial_binding_sha256': sha256(args.output / 'initial_binding.json'),
        'root_decision_sha256': args.decision_sha, 'template_sha256': TEMPLATE_SHA,
        'output_bindings': {name: sha256(args.output / name) for name in [
            'safs/inputs.json', 'safs/protocol.json', 'window_comparison/protocol.json', 'window_comparison/comparison.json']},
        'draft_sha256': sha256(draft_path), 'submitted': False,
        'strict_qualification': False, 'scientifically_ready': False,
        'development_selection_disclosed': True,
        'known_development_decline_fraction_relative_error': .4120711578087457,
        'scientific_qualification': 'Best available development-selected fallback; strict recovery and seed-stability gates failed.'}
    write_json(args.output / 'execution_binding.json', result)
    print(json.dumps({'draft_created': True, 'strict_qualification': False,
                      'failed_gates': [k for k, passed in gates.items() if not passed]}))


if __name__ == '__main__':
    main()
