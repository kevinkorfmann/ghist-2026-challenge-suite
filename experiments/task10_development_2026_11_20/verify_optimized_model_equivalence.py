"""Compare frozen-model probabilities on retained training features only."""
import json
import os
from pathlib import Path

import joblib
import numpy as np

from scripts.develop_task10_detected_ibd import save_json, sha256, verify_inputs


def main():
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('Use Slurm')
    root = Path('/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task10-development-20261120')
    evaluation = root / 'detected_panel_v2/model_evaluation'
    optimized = root / 'optimized_scaling_v1'
    output = root / 'optimized_model_equivalence_v1'
    output.mkdir(exist_ok=False)
    freeze = json.loads((evaluation / 'model_calibration_freeze.json').read_text())
    equivalent = json.loads((optimized / 'equivalence.json').read_text())
    verify_inputs(equivalent['input_sha256'])
    checks, inputs = [], dict(equivalent['input_sha256'])
    inputs[str(evaluation / 'model_calibration_freeze.json')] = sha256(evaluation / 'model_calibration_freeze.json')
    for variant in ['context', 'pair_only']:
        path = evaluation / f'{variant}_model.joblib'
        if sha256(path) != freeze['model_hashes'][variant]:
            raise ValueError('Final frozen model changed')
        inputs[str(path)] = sha256(path)
        bundle = joblib.load(path)
        model = bundle['model']
        # Serial accumulation avoids nondeterministic thread reduction roundoff.
        model.n_jobs = 1
        for entry in equivalent['checks']:
            index = entry['cohort_index']
            original = json.loads((root / f'detected_panel_v2/cohort_{index:03d}/features.json').read_text())
            fast_path = optimized / f'training_{index:03d}_features.json'
            faster = json.loads(fast_path.read_text())
            inputs[str(fast_path)] = sha256(fast_path)
            matrices = [np.asarray([[row[name] for name in bundle['feature_names']] for row in rows], dtype=np.float32)
                for rows in [original, faster]]
            if not np.array_equal(*matrices):
                raise ValueError('Actual float32 model inputs differ')
            predictions = [model.predict_proba(matrix) for matrix in matrices]
            if not np.array_equal(*predictions):
                raise ValueError('Model probabilities differ')
            checks.append({'variant': variant, 'training_cohort_index': index, 'pairs': len(original),
                'float32_model_inputs_identical': True, 'serial_probabilities_identical': True})
    verify_inputs(inputs)
    save_json(output / 'result.json', {'state': 'MODEL_EQUIVALENCE_VERIFIED', 'checks': checks,
        'input_sha256': inputs, 'source_sha256': sha256(Path(__file__)),
        'held_out_outcomes_read': False, 'biological_qualification': False})


if __name__ == '__main__':
    main()
