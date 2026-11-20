"""A single model-implied-distance tree diagnostic; no fitting or true topology."""
from __future__ import annotations

import argparse
import itertools
import json
import math
import os
from pathlib import Path

from calibrate_phase import metrics
from ghist_relatedness.labels import HUMAN_LABELS

from scripts.develop_task10_detected_ibd import global_tree_predictions, save_json, sha256, verify_inputs


def model_distances(predictions):
    """Encode the forest's mean meiosis count in the frozen decoder's distance API.

    The existing decoder transforms r to d=1-log2(r). Encoding r=2**(1-E[d])
    therefore changes only the inferred distance input, retaining its NJ, integer
    edge, valid-root, score, tie and fallback policies exactly. These synthetic
    coefficients are not measured IBD or kinship, and imply no calibrated root
    posterior.
    """
    meioses = {label: up + down for (up, down, _), label in HUMAN_LABELS.items()}
    samples = sorted({str(row[key]) for row in predictions for key in ['indID1', 'indID2']})
    pairs, features, distances = set(), [], []
    for row in predictions:
        pair = tuple(sorted([row['indID1'], row['indID2']]))
        if pair in pairs or pair[0] == pair[1]:
            raise ValueError('Duplicate or self pair')
        pairs.add(pair)
        probabilities = row['probabilities']
        if set(probabilities) != set(meioses) or any(not math.isfinite(p) or p < 0 for p in probabilities.values()):
            raise ValueError('Invalid complete-label probabilities')
        if not math.isclose(sum(probabilities.values()), 1., rel_tol=1e-10, abs_tol=1e-10):
            raise ValueError('Probabilities must sum to one')
        distance = math.fsum(probabilities[label] * value for label, value in meioses.items())
        features.append({'indID1': pair[0], 'indID2': pair[1], 'relatedness_ge_10': 2 ** (1 - distance)})
        distances.append({'indID1': pair[0], 'indID2': pair[1], 'model_mean_meioses': distance})
    if pairs != set(itertools.combinations(samples, 2)):
        raise ValueError('Incomplete roster pairs')
    return features, distances


def decode(predictions):
    features, distances = model_distances(predictions)
    rows, tree = global_tree_predictions(features, predictions, minimum_cm=10)
    tree.update(distance_source='Frozen forest mean meioses, not measured IBD or kinship',
        threshold_cm=None, model_implied_distances=distances, uncertainty_calibrated=False)
    return rows, tree


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', type=Path, required=True)
    parser.add_argument('--protocol', type=Path, required=True)
    parser.add_argument('--panel-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):
        parser.error('Use Slurm')
    protocol = json.loads(args.protocol.read_text())
    if sha256(Path(__file__)) != protocol['decoder_source_sha256']:
        raise ValueError('Decoder source changed')
    if sha256(args.source_root / 'source_manifest.json') != protocol['candidate_source_manifest_sha256']:
        raise ValueError('Frozen source changed')
    verify_inputs(json.loads((args.source_root / 'source_manifest.json').read_text()))
    evaluation = args.panel_root / 'model_evaluation'
    freeze = evaluation / 'model_calibration_freeze.json'
    result_path = evaluation / 'result.json'
    result = json.loads(result_path.read_text())
    if result['advance_to_confirmation']:
        raise ValueError('This follow-up is only for failed v2 development')
    bound = {str(path): sha256(path) for path in [args.protocol, Path(__file__), freeze, result_path]}
    args.output.mkdir(parents=True, exist_ok=False)
    scores, diagnostics = {}, {}
    for variant in protocol['variants']:
        source_path = evaluation / f'{variant}_development_scorer_input.json'
        source = json.loads(source_path.read_text())
        bound[str(source_path)] = sha256(source_path)
        cohorts, trees = [], []
        for cohort in source['cohorts']:
            predictions, tree = decode(cohort['predictions'])
            # Truth is attached only after inference and never passed to decode.
            cohorts.append({'cohort_id': cohort['cohort_id'], 'truth': cohort['truth'], 'predictions': predictions})
            trees.append({'cohort_id': cohort['cohort_id'], **tree})
        save_json(args.output / f'{variant}_scorer_input.json', {'labels': source['labels'], 'cohorts': cohorts})
        save_json(args.output / f'{variant}_predictions.json', [{**row, 'cohort_id': cohort['cohort_id']}
            for cohort in cohorts for row in cohort['predictions']])
        save_json(args.output / f'{variant}_trees.json', trees)
        scores[variant] = metrics(cohorts, source['labels'])
        diagnostics[variant] = {'valid_root_counts': {row['cohort_id']: row['valid_roots'] for row in trees},
            'fallback_cohorts': [row['cohort_id'] for row in trees if not row['valid_roots']]}
    verify_inputs(bound)
    save_json(args.output / 'result.json', {'stage': 'POST_V2_FAILURE_DEVELOPMENT_DIAGNOSTIC',
        'protocol_sha256': sha256(args.protocol), 'scores': scores, 'diagnostics': diagnostics,
        'input_sha256': bound, 'original_v2_policy_changed': False, 'uncertainty_calibrated': False,
        'submission_ready': False, 'advance_to_confirmation': False})


if __name__ == '__main__':
    main()
