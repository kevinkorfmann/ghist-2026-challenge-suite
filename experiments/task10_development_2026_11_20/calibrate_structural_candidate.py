"""Apply one predeclared global correctness rule to retained calibration families."""
from __future__ import annotations

import json
import os
from pathlib import Path

import joblib
from structural_confidence import attach_global_confidence, global_correctness, reliability
from structural_pair_candidate import infer, pair_predictions

from scripts.develop_task10_detected_ibd import save_json, sha256, verify_inputs


def main():
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('Use Slurm')
    root = Path('/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task10-development-20261120')
    source = Path(__file__).resolve().parent
    manifest = json.loads((source / 'source_manifest.json').read_text())
    verify_inputs(manifest)
    protocol = json.loads((source / 'structural_qualification_protocol_v3.json').read_text())
    evaluation = root / 'detected_panel_v2/model_evaluation'
    freeze_path = evaluation / 'model_calibration_freeze.json'
    freeze = json.loads(freeze_path.read_text())
    verify_inputs(freeze['hashes'])
    model_path = evaluation / 'pair_only_model.joblib'
    if sha256(model_path) != freeze['model_hashes']['pair_only']:
        raise ValueError('Frozen pair model changed')
    bundle = joblib.load(model_path)
    labels = list(map(str, bundle['model'].classes_))
    output = root / 'structural_calibration_v3'
    output.mkdir(exist_ok=False)
    inputs = {str(path): sha256(path) for path in [model_path, freeze_path]}
    inferred = []
    for index in protocol['confidence_calibration']['cohort_indices']:
        cohort = root / f'detected_panel_v2/cohort_{index:03d}'
        path = cohort / 'features.json'
        if sha256(path) != freeze['hashes'][str(path)]:
            raise ValueError('Retained calibration features changed')
        inputs[str(path)] = sha256(path)
        predictions = pair_predictions(json.loads(path.read_text()), bundle,
            freeze['temperatures']['pair_only']['temperature'])
        rows, tree = infer(predictions)
        save_json(output / f'cohort_{index:03d}_inference.json', {'predictions': rows, 'tree': tree})
        inferred.append({'cohort_id': f'detected_panel_{index:03d}', 'predictions': rows})
    # All inference is complete and retained before reading calibration outcomes.
    for index, cohort in zip(protocol['confidence_calibration']['cohort_indices'], inferred, strict=True):
        path = root / f'detected_panel_v2/cohort_{index:03d}/truth.json'
        if sha256(path) != freeze['hashes'][str(path)]:
            raise ValueError('Retained calibration truth changed')
        inputs[str(path)] = sha256(path)
        cohort['truth'] = json.loads(path.read_text())
    calibration = global_correctness(inferred, labels)
    if calibration['pairs'] != protocol['confidence_calibration']['expected_pairs']:
        raise ValueError('Incomplete predeclared calibration coverage')
    calibration.update(protocol_sha256=sha256(source / 'structural_qualification_protocol_v3.json'),
        source_manifest_sha256=sha256(source / 'source_manifest.json'), input_sha256=inputs,
        qualification_outcomes_observed=False)
    save_json(output / 'confidence_freeze.json', calibration)
    cohorts = [{**cohort, 'predictions': attach_global_confidence(cohort['predictions'], labels, calibration)}
        for cohort in inferred]
    save_json(output / 'scorer_input.json', {'labels': labels, 'cohorts': cohorts})
    save_json(output / 'reliability.json', reliability(cohorts, labels))
    verify_inputs(inputs)
    verify_inputs(manifest)


if __name__ == '__main__':
    main()
