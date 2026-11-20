"""Fixed development-only scalar calibration and advancement checks."""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np


def digest(path):
    value = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def key(row):
    return tuple(sorted((row['indID1'], row['indID2'])))


def metrics(cohorts, labels):
    counts, correct = Counter(), Counter()
    by_cohort = {}
    for cohort in cohorts:
        truth = {key(row): row['relation'] for row in cohort['truth']}
        predictions = {key(row): row['relation'] for row in cohort['predictions']}
        if len(truth) != len(cohort['truth']) or len(predictions) != len(cohort['predictions']):
            raise ValueError('Duplicate pair')
        if set(truth) != set(predictions):
            raise ValueError('Missing or extra pair')
        if any(a == b for a, b in truth):
            raise ValueError('Self pair')
        if (set(truth.values()) | set(predictions.values())) - set(labels):
            raise ValueError('Unsupported label')
        cohort_correct = 0
        for pair, label in truth.items():
            counts[label] += 1
            correct[label] += predictions[pair] == label
            cohort_correct += predictions[pair] == label
        by_cohort[cohort['cohort_id']] = cohort_correct / len(truth)
    recall = {label: correct[label] / counts[label] if counts[label] else None for label in labels}
    return {'pairs': sum(counts.values()), 'correct': sum(correct.values()),
        'exact_accuracy': sum(correct.values()) / sum(counts.values()),
        'macro_recall': float(np.mean([recall[label] or 0 for label in labels])),
        'minimum_label_recall': min(recall[label] or 0 for label in labels),
        'recall': recall, 'counts': dict(counts), 'cohort_accuracy': by_cohort,
        'complete_label_coverage': all(counts[label] for label in labels)}


def temperature(probabilities, value):
    logs = np.log(np.maximum(probabilities, 1e-12)) / value
    logs -= logs.max(axis=1, keepdims=True)
    result = np.exp(logs)
    return result / result.sum(axis=1, keepdims=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--experiment-root', type=Path, required=True)
    parser.add_argument('--out-dir', type=Path, required=True)
    args = parser.parse_args()
    root, output = args.experiment_root.resolve(), args.out_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    protocol = json.loads((root / 'phase_protocol.json').read_text())
    training = json.loads((root / 'phase_training/result.json').read_text())
    labels = training['labels']
    inputs = {str(Path(__file__).resolve()): digest(Path(__file__)),
        str(root / 'phase_protocol.json'): digest(root / 'phase_protocol.json')}
    cohorts = {variant: [] for variant in ['phase', 'phase_free', 'baseline']}
    for variant in ['phase', 'phase_free']:
        model = root / f'phase_training/{variant}_model.joblib'
        if digest(model) != training['model_hashes'][variant]:
            raise ValueError('Trained model bytes changed')
        inputs[str(model)] = digest(model)
    for index in range(4):
        directory = root / f'phase_development_cohort_{index}'
        truth_path = directory / 'truth.json'
        truth = json.loads(truth_path.read_text())
        inputs[str(truth_path)] = digest(truth_path)
        for variant, variant_cohorts in cohorts.items():
            if variant == 'baseline':
                path = directory / 'representation_comparison/raw_haplotype.json'
            else:
                path = root / f'phase_predictions/cohort_{index}/{variant}/predictions.json'
                run_protocol = json.loads((path.parent / 'protocol.json').read_text())
                for bound, expected in run_protocol['hashes'].items():
                    if digest(Path(bound)) != expected:
                        raise ValueError(f'Prediction input changed: {bound}')
                if run_protocol['hashes'][str(truth_path)] != inputs[str(truth_path)]:
                    raise ValueError('Scoring truth differs from prediction-bound truth')
            predictions = json.loads(path.read_text())
            inputs[str(path)] = digest(path)
            variant_cohorts.append({'cohort_id': f'phase_development_{index}',
                'truth': truth, 'predictions': predictions})
    scores = {variant: metrics(rows, labels) for variant, rows in cohorts.items()}
    values, truth_labels = [], []
    for cohort in cohorts['phase']:
        expected = {key(row): row['relation'] for row in cohort['truth']}
        for row in cohort['predictions']:
            probabilities = row['probabilities']
            if set(probabilities) != set(labels):
                raise ValueError('Incomplete probability domain')
            values.append([probabilities[label] for label in labels])
            truth_labels.append(labels.index(expected[key(row)]))
    values, truth_labels = np.asarray(values), np.asarray(truth_labels)
    if not np.isfinite(values).all() or np.any(values < 0) or not np.allclose(values.sum(axis=1), 1):
        raise ValueError('Invalid prediction probabilities')
    counts = np.bincount(truth_labels, minlength=len(labels))
    if np.any(counts == 0):
        raise ValueError('Cannot calibrate across missing labels')
    weights = 1 / counts[truth_labels]
    losses = []
    for value in protocol['scalar_temperature_grid']:
        adjusted = temperature(values, value)
        loss = -float(np.sum(weights * np.log(np.maximum(adjusted[np.arange(len(values)), truth_labels], 1e-12))) / weights.sum())
        losses.append({'temperature': value, 'class_balanced_log_loss': loss})
    selected = min(losses, key=lambda row: (row['class_balanced_log_loss'], row['temperature']))
    score, baseline = scores['phase'], scores['baseline']
    gates = {'all_labels': score['complete_label_coverage'], 'exact_accuracy': score['exact_accuracy'] >= .8,
        'macro_recall': score['macro_recall'] >= .7, 'minimum_label_recall': score['minimum_label_recall'] >= .5,
        'gain_over_baseline': score['exact_accuracy'] - baseline['exact_accuracy'] >= .1,
        'no_cohort_regression': all(value >= baseline['cohort_accuracy'][cohort] for cohort, value in score['cohort_accuracy'].items())}
    for path, expected in inputs.items():
        if digest(Path(path)) != expected:
            raise ValueError('Bound input changed during calibration')
    for variant, rows in cohorts.items():
        save(output / f'{variant}_scorer_input.json', {'labels': labels, 'cohorts': rows})
    save(output / 'result.json', {'stage': 'DEVELOPMENT', 'scores': scores, 'gates': gates,
        'advance_to_confirmation': all(gates.values()), 'temperature_selection': selected,
        'temperature_grid_results': losses, 'inputs_sha256': inputs,
        'calibration_interpretation': 'Selected on these development cohorts; not independent calibration validation',
        'submission_ready': False, 'submission_count': 0})
    print(json.dumps({'scores': scores, 'gates': gates, 'temperature': selected}, indent=2))


if __name__ == '__main__':
    main()
