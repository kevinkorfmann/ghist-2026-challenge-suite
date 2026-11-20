"""Inference-only adapter for the frozen pair forest and integer-edge decoder."""
from __future__ import annotations

import itertools

import numpy as np
from calibrate_phase import temperature
from forest_distance_decoder import decode
from integer_edge_decoder import optimize

from scripts.develop_task10_detected_ibd import phase_pair_features


def pair_features(samples, segments, genetic_map, genome_cm):
    if len(set(samples)) != len(samples):
        raise ValueError('Duplicate sample')
    return [{'indID1': a, 'indID2': b,
        **phase_pair_features(segments.get((a, b), {}), genetic_map, genome_cm)}
        for a, b in itertools.combinations(sorted(samples), 2)]


def pair_predictions(features, bundle, calibration_temperature):
    if any(name.startswith('context_') for name in bundle['feature_names']):
        raise ValueError('Expected frozen pair-only model')
    matrix = np.asarray([[row[name] for name in bundle['feature_names']] for row in features], dtype=np.float32)
    # Serial forest reduction fixes accumulation order without changing any tree.
    bundle['model'].n_jobs = 1
    values = temperature(bundle['model'].predict_proba(matrix), calibration_temperature)
    return [{key: row[key] for key in ['indID1', 'indID2']} |
        {'relation': str(bundle['model'].classes_[probabilities.argmax()]),
         'probabilities': dict(zip(map(str, bundle['model'].classes_), map(float, probabilities), strict=True))}
        for row, probabilities in zip(features, values, strict=True)]


def infer(predictions):
    _, initial = decode(predictions)
    if not initial['valid_roots']:
        raise ValueError('Initial NJ tree has no accepted root; candidate cannot provide structural predictions')
    rows, tree = optimize(predictions, initial['integer_edges'], max_moves=20, max_feasibility_checks=200)
    return rows, {'initial_tree': initial, 'optimized_tree': tree,
        'uncertainty_calibrated': False, 'biological_qualification': False}
