"""Global empirical correctness calibration; no individual or root posterior."""
from __future__ import annotations

import math


def global_correctness(cohorts, labels):
    correct, count = 0, 0
    observed = set()
    for cohort in cohorts:
        def keys(rows):
            return {tuple(sorted([row['indID1'], row['indID2']])): row['relation'] for row in rows}
        truth, predictions = keys(cohort['truth']), keys(cohort['predictions'])
        if len(truth) != len(cohort['truth']) or len(predictions) != len(cohort['predictions']) or set(truth) != set(predictions):
            raise ValueError('Incomplete or duplicated calibration pairs')
        if (set(truth.values()) | set(predictions.values())) - set(labels):
            raise ValueError('Unsupported calibration label')
        observed.update(truth.values())
        correct += sum(predictions[pair] == label for pair, label in truth.items())
        count += len(truth)
    if observed != set(labels) or not count:
        raise ValueError('Calibration must cover every label')
    q = (correct + 1) / (count + 2)
    if q <= 1 / len(labels):
        raise ValueError('Global correctness cannot support the selected top label')
    return {'correct': correct, 'pairs': count, 'confidence': q,
        'rule': '(correct + 1) / (pairs + 2), all fixed calibration pairs; no tuning',
        'individual_resolution': 0, 'root_posterior': False,
        'remaining_label_rule': 'Uniform (1-confidence)/(label_count-1); no discrimination among alternatives'}


def attach_global_confidence(predictions, labels, calibration):
    q = calibration['confidence']
    if len(set(labels)) != len(labels) or not math.isfinite(q) or not 1 / len(labels) < q < 1:
        raise ValueError('Invalid global correctness confidence')
    if any(row['relation'] not in labels for row in predictions):
        raise ValueError('Unsupported predicted label')
    return [{**row, 'probabilities': {label: q if label == row['relation'] else (1 - q) / (len(labels) - 1)
        for label in labels}} for row in predictions]


def reliability(cohorts, labels):
    """Top-label ECE and transparent heterogeneity by true/predicted label/cohort."""
    groups = {'all': [], **{f'truth:{label}': [] for label in labels},
        **{f'predicted:{label}': [] for label in labels}}
    for cohort in cohorts:
        name = f'cohort:{cohort["cohort_id"]}'
        groups[name] = []
        truth = {tuple(sorted([row['indID1'], row['indID2']])): row['relation'] for row in cohort['truth']}
        if len(truth) != len(cohort['truth']) or len(cohort['predictions']) != len(truth):
            raise ValueError('Incomplete reliability pair set')
        seen = set()
        for row in cohort['predictions']:
            pair = tuple(sorted([row['indID1'], row['indID2']]))
            if pair in seen or pair not in truth:
                raise ValueError('Extra or duplicate reliability pair')
            seen.add(pair)
            values = row['probabilities']
            if set(values) != set(labels) or any(not math.isfinite(p) or p < 0 for p in values.values()):
                raise ValueError('Invalid complete probability vector')
            if not math.isclose(sum(values.values()), 1., abs_tol=1e-10) or values[row['relation']] != max(values.values()):
                raise ValueError('Top-label confidence contract failed')
            value = (float(row['relation'] == truth[pair]), values[row['relation']])
            for key in ['all', name, f'truth:{truth[pair]}', f'predicted:{row["relation"]}']:
                groups[key].append(value)
    result = {}
    for name, rows in groups.items():
        if not rows:
            result[name] = {'pairs': 0, 'accuracy': None, 'mean_confidence': None, 'ece_10bins': None}
            continue
        bins = [[] for _ in range(10)]
        for actual, confidence in rows:
            bins[min(int(confidence * 10), 9)].append((actual, confidence))
        ece = sum(abs(sum(actual - confidence for actual, confidence in values)) for values in bins) / len(rows)
        result[name] = {'pairs': len(rows), 'accuracy': sum(a for a, _ in rows) / len(rows),
            'mean_confidence': sum(p for _, p in rows) / len(rows), 'ece_10bins': ece}
    return result
