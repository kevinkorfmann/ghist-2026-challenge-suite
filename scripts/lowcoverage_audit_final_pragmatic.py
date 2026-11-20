"""Independent arithmetic/schema audit; never estimates, downloads, or submits.

Only the root agent should invoke the CLI against actual challenge outputs.
This helper was developed using archived synthetic spectra and fixtures only.
"""
import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import re

import numpy as np
import yaml


FREEZE = '7a98a0e446cf5b5221fef48bace8743b5a604f7cd807fdff0c26df3df001b1bb'
TEMPLATES = {
    'split_isolation': ('split_isolation_low_coverage.submission.yaml', '46bddd26dd3ef9fac2d707c549214b5b8f9d2c427b9d0a38cc7ecc9a4ed2d0cb'),
    'bottleneck': ('bottleneck.submission.yaml', 'd2d4ad595d4de0466fc4291d3f6550054da98676e42c662b6b57aaeabf5a9936'),
}


def digest(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def close(actual, expected, name):
    if not math.isfinite(float(actual)) or not math.isclose(actual, expected, rel_tol=1e-10, abs_tol=1e-8):
        raise ValueError(f'{name} mismatch: {actual} != {expected}')


def spectrum(path, shape):
    text = Path(path).read_text()
    match = re.match(r'^#SHAPE=<([0-9/, ]+)>\s*\n', text)
    if not match or [int(x) for x in re.split(r'[,/ ]+', match[1])] != shape:
        raise ValueError('spectrum shape mismatch')
    values = np.array(text[match.end():].split(), dtype=float)
    if values.size != math.prod(shape) or not np.all(np.isfinite(values)) or np.any(values < 0):
        raise ValueError('spectrum must contain the complete nonnegative finite count array')
    return values


def audit_comparison(directory, task='split_isolation'):
    directory = Path(directory)
    split = task == 'split_isolation'
    shape, cap, tolerance, window = ([21, 21], 500, 1e-6, 500) if split else ([41], 50, 1e-4, 100)
    protocol = read(directory / 'protocol.json')
    if (protocol['sample_sizes'] != [n - 1 for n in shape] or protocol['blocks'] != 500
            or protocol['window_size'] != window or protocol['max_epochs'] != cap
            or protocol['tolerance'] != tolerance):
        raise ValueError('unexpected frozen estimator policy')
    recorded = read(directory / 'comparison.json')
    arrays, targets, runs = {}, {}, {}
    bounds = [(0.05, 20), (0.05, 20), (0.0001, 10)] if split else [(0.001, 0.9999), (0.0001, 10)]
    for label in ['seed1', 'seed2', 'continued']:
        values = spectrum(directory / f'{label}.folded.sfs', shape)
        unfolded = spectrum(directory / f'{label}.unfolded.sfs', shape).reshape(shape)
        expected = np.zeros(shape)
        half = sum(n - 1 for n in shape) / 2
        for index in np.ndindex(*shape):
            complement = tuple(n - 1 - i for n, i in zip(shape, index))
            if sum(index) < half:
                expected[index] = unfolded[index] + unfolded[complement]
            elif sum(index) == half:
                expected[index] = (unfolded[index] + unfolded[complement]) / 2
        if not np.allclose(values, expected.ravel(), rtol=0, atol=1.1e-6):
            raise ValueError(f'{label}: global folding differs from raw unfolded counts')
        length = float(values.sum())
        if not 0 < length <= 10_000_001:
            raise ValueError('invalid full-width exposure')
        estimate = read(directory / f'{label}.estimate.json')
        records = estimate['restarts']
        if Counter(row['method'] for row in records) != {'L-BFGS-B': 8, 'Powell': 8}:
            raise ValueError('sixteen optimizer records required')
        for row in records:
            parameters = np.asarray(row['parameters'])
            if parameters.shape != (len(bounds),) or not np.all(np.isfinite(parameters)) or np.any(parameters <= 0):
                raise ValueError('invalid optimizer parameters')
            boundary = any(abs(math.log(p / lo)) < 1e-3 or abs(math.log(p / hi)) < 1e-3
                           for p, (lo, hi) in zip(parameters, bounds))
            if boundary != row['at_boundary']:
                raise ValueError('optimizer boundary flag mismatch')
        finite = [row for row in records if math.isfinite(row['negative_log_likelihood'])
                  and row['negative_log_likelihood'] < 1e90]
        best = min(finite, key=lambda row: row['negative_log_likelihood'])
        if best != estimate['best_fit']:
            raise ValueError('reported best fit is not global finite best')
        agreeing = [row for row in finite if row['success'] and not row['at_boundary']
                    and abs(row['negative_log_likelihood'] - best['negative_log_likelihood']) < .01
                    and max(abs(math.log(p / q)) for p, q in zip(row['parameters'], best['parameters'])) < .01]
        optimizer_ok = bool(best['success'] and not best['at_boundary']
                            and {row['method'] for row in agreeing} == {'L-BFGS-B', 'Powell'})
        if optimizer_ok != estimate['optimizer_agreement']['passed']:
            raise ValueError('optimizer agreement mismatch')
        theta = best['theta']
        if not math.isfinite(theta) or theta <= 0:
            raise ValueError('invalid theta')
        ancestor = theta / (4e-7 * length)
        pars = best['parameters']
        physical = ({'pop1_population_size': ancestor * pars[0],
                     'pop2_population_size': ancestor * pars[1],
                     'generations_since_split': 2 * ancestor * pars[2]} if split else
                    {'post_decline_fraction': pars[0], 'generations': 2 * ancestor * pars[1]})
        close(estimate['ancestral_population_size'], ancestor, 'ancestral size')
        close(estimate['callable_length'], length, 'callable length')
        close(estimate['callable_fraction'], length / 1e7, 'callable fraction')
        for name, value in physical.items():
            close(estimate['targets'][name], value, name)
            close(recorded['estimates'][label][name], value, name)
        log = (directory / f'{label}.log').read_text()
        epochs = [int(value) for value in re.findall(r'Finished epoch ([0-9]+)', log)]
        stops = re.findall(r'Δ=([-+0-9.e]+)\s*([≤>])\s*([-+0-9.e]+)', log)
        if not epochs or not stops:
            raise ValueError('missing actual winsfs epoch/stopping diagnostics')
        delta, relation, declared = stops[-1]
        close(float(declared), tolerance, 'stopping tolerance')
        stopped = epochs[-1] < cap and relation == '≤' and abs(float(delta)) <= tolerance
        if estimate['winsfs_epochs'] != epochs or estimate['winsfs_early_stop_observed'] != (epochs[-1] < cap):
            raise ValueError('logged epochs disagree with estimate')
        arrays[label], targets[label] = values, physical
        runs[label] = {'targets': physical, 'callable_fraction': length / 1e7,
                       'optimizer_agreement': optimizer_ok, 'epochs': epochs[-1],
                       'last_delta': float(delta), 'stopping_passed': stopped}
    comparisons = {}
    for label in ['seed2', 'continued']:
        a, b = arrays['seed1'], arrays[label]
        if a[1:].sum() <= 0:
            raise ValueError('no segregating information')
        row = {'frequency_L1': float(np.abs(a / a.sum() - b / b.sum()).sum()),
               'segregating_relative_L1': float(np.abs(a[1:] - b[1:]).sum() / a[1:].sum()),
               'target_relative_change': {k: targets[label][k] / v - 1 for k, v in targets['seed1'].items()}}
        for key in ['frequency_L1', 'segregating_relative_L1']:
            close(recorded['comparisons'][label][key], row[key], key)
        for key, value in row['target_relative_change'].items():
            close(recorded['comparisons'][label]['target_relative_change'][key], value, key)
        comparisons[label] = row
    gates = {'exposure': all(row['callable_fraction'] >= .99 for row in runs.values()),
             'observed_stopping': all(row['stopping_passed'] for row in runs.values()),
             'two_optimizer_agreement': all(row['optimizer_agreement'] for row in runs.values()),
             'target_stability': all(abs(v) < .01 for row in comparisons.values() for v in row['target_relative_change'].values()),
             'frequency_L1': all(row['frequency_L1'] < .001 for row in comparisons.values()),
             'segregating_relative_L1': all(row['segregating_relative_L1'] < .05 for row in comparisons.values())}
    return {'task': task, 'runs': runs, 'comparisons': comparisons, 'strict_numerical_gates': gates,
            'practical_5percent_target_stability': all(abs(v) <= .05 for row in comparisons.values()
                                                     for v in row['target_relative_change'].values()),
            'strict_qualification': False, 'truth_or_recovery_checked_on_final': False}


def audit_draft(path, template, targets, task):
    if digest(template) != TEMPLATES[task][1]:
        raise ValueError('official template hash mismatch')
    actual, expected = yaml.safe_load(Path(path).read_text()), yaml.safe_load(Path(template).read_text())
    if set(actual) != set(expected) or set(actual['parameters']) != set(targets):
        raise ValueError('official YAML schema mismatch')
    for key in expected:
        if key != 'parameters' and actual[key] != expected[key]:
            raise ValueError('official template metadata changed')
    for name, value in actual['parameters'].items():
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
            raise ValueError('draft parameter is not a positive finite number')
        if name == 'post_decline_fraction' and value >= 1:
            raise ValueError('decline fraction must be below one')
        close(value, targets[name], name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--frozen-root', type=Path, required=True)
    parser.add_argument('--task', choices=list(TEMPLATES), default='split_isolation')
    parser.add_argument('--root-decision', type=Path, required=True)
    parser.add_argument('--decision-sha', required=True)
    args = parser.parse_args()
    output, root = args.output, args.frozen_root
    result = audit_comparison(output / 'window_comparison', args.task)
    initial, execution = read(output / 'initial_binding.json'), read(output / 'execution_binding.json')
    if digest(args.root_decision) != args.decision_sha:
        raise ValueError('root decision hash mismatch')
    decision = read(args.root_decision)
    primary = 'seed1' if args.task == 'split_isolation' else 'continued'
    if not (initial['freeze_sha256'] == execution['freeze_sha256'] == FREEZE
            and initial['root_decision_sha256'] == execution['root_decision_sha256'] == args.decision_sha
            and initial['wrapper_sha256'] == execution['wrapper_sha256'] == decision['wrapper_sha256']
            and initial['primary_estimate'] == decision['primary_estimate'] == primary
            and decision['authorize_final_inference'] is True
            and execution['source_runtime_input_pre_post_verified'] is True
            and execution['strict_qualification'] is False and execution['scientifically_ready'] is False
            and initial['input_binding']['audit'] == read(output / 'safs/inputs.json')):
        raise ValueError('execution authorization/input bindings mismatch')
    if digest(output / 'initial_binding.json') != execution['initial_binding_sha256']:
        raise ValueError('initial binding hash mismatch')
    for name, expected in execution['output_bindings'].items():
        if digest(output / name) != expected:
            raise ValueError(f'output hash mismatch: {name}')
    for name, value in execution['primary_estimate'].items():
        close(value, result['runs'][primary]['targets'][name], name)
    draft = output / ('task14.pragmatic.submission.draft.yaml' if args.task == 'split_isolation'
                      else 'task13.pragmatic.submission.draft.yaml')
    audit_draft(draft, root / 'official_templates' / TEMPLATES[args.task][0], result['runs'][primary]['targets'], args.task)
    if digest(draft) != execution['draft_sha256']:
        raise ValueError('draft hash mismatch')
    result.update({'primary_label': primary, 'draft_sha256': digest(draft),
                   'execution_binding_sha256': digest(output / 'execution_binding.json'),
                   'root_decision_sha256': args.decision_sha,
                   'reported_numerical_gates': execution['numerical_gates'],
                   'numeric_schema_binding_audit_passed': True,
                   'raw_input_and_runtime_rehash': 'Root must separately retain pre/post hash audit; this helper does not repeat large input/runtime reads.'})
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))


if __name__ == '__main__':
    main()
