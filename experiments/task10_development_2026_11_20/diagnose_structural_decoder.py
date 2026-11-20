"""One fixed development-only NJ10 decoder of frozen context or pair probabilities."""
import argparse
import json
import os
from pathlib import Path

from calibrate_phase import metrics

from scripts.develop_task10_detected_ibd import global_tree_predictions, save_json, sha256, verify_inputs


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
        raise ValueError('Frozen input source changed')
    verify_inputs(json.loads((args.source_root / 'source_manifest.json').read_text()))
    evaluation = args.panel_root / 'model_evaluation'
    # The original fit/calibration/evaluation must have finished and cannot be changed.
    result_path = evaluation / 'result.json'
    result = json.loads(result_path.read_text())
    if result['stage'] != 'DEVELOPMENT':
        raise ValueError('Unexpected input stage')
    freeze = evaluation / 'model_calibration_freeze.json'
    bound = {str(result_path): sha256(result_path), str(freeze): sha256(freeze),
        str(args.protocol): sha256(args.protocol), str(Path(__file__)): sha256(Path(__file__))}
    args.output.mkdir(parents=True, exist_ok=False)
    if result['advance_to_confirmation']:
        verify_inputs(bound)
        save_json(args.output / 'result.json', {'state': 'SKIPPED_PRIMARY_V2_PASSED',
            'input_sha256': bound, 'submission_ready': False})
        return
    scores, diagnostics = {}, {}
    for variant in protocol['variants']:
        prediction_path = evaluation / f'{variant}_development_scorer_input.json'
        source = json.loads(prediction_path.read_text())
        bound[str(prediction_path)] = sha256(prediction_path)
        cohorts, trees = [], []
        for cohort in source['cohorts']:
            index = int(cohort['cohort_id'].rsplit('_', 1)[1])
            directory = args.panel_root / f'cohort_{index:03d}'
            record = json.loads((directory / 'result.json').read_text())
            if record['item']['stage'] != 'development':
                raise ValueError('Unexpected cohort split')
            feature_path = directory / 'features.json'
            bound[str(feature_path)] = record['output_sha256'][str(feature_path)]
            verify_inputs({str(feature_path): bound[str(feature_path)]})
            rows = json.loads(feature_path.read_text())
            features = [{key: row[key] for key in ['indID1', 'indID2']} |
                {'relatedness_ge_10': row['state_relatedness_ge_10']} for row in rows]
            predictions, diagnostic = global_tree_predictions(features, cohort['predictions'], minimum_cm=10)
            # Root score is composite; no probability or certainty is exported.
            cohorts.append({'cohort_id': cohort['cohort_id'], 'truth': cohort['truth'], 'predictions': predictions})
            trees.append({'cohort_id': cohort['cohort_id'], **diagnostic})
        save_json(args.output / f'{variant}_scorer_input.json', {'labels': source['labels'], 'cohorts': cohorts})
        save_json(args.output / f'{variant}_predictions.json', [{**row, 'cohort_id': cohort['cohort_id']}
            for cohort in cohorts for row in cohort['predictions']])
        save_json(args.output / f'{variant}_trees.json', trees)
        scores[variant] = metrics(cohorts, source['labels'])
        diagnostics[variant] = {'valid_root_counts': {row['cohort_id']: row['valid_roots'] for row in trees},
            'fallback_cohorts': [row['cohort_id'] for row in trees if not row['valid_roots']]}
    verify_inputs(bound)
    save_json(args.output / 'result.json', {'stage': 'SUPPLEMENTARY_DEVELOPMENT_ONLY',
        'protocol_sha256': sha256(args.protocol), 'scores': scores, 'diagnostics': diagnostics,
        'input_sha256': bound, 'original_v2_gates': result['gates'],
        'original_v2_policy_changed': False, 'uncertainty_calibrated': False,
        'submission_ready': False, 'advance_to_confirmation': False})


if __name__ == '__main__':
    main()
