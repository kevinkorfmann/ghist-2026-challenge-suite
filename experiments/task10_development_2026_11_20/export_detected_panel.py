"""Export frozen detected-panel development predictions with raw provenance."""
from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path


def main():
    root = Path(__file__).resolve().parent
    source = Path('/tmp/task10-evidence-20261120/detected_panel_v2')
    destination = root / 'detected_panel_evidence'
    if not (source / 'model_evaluation/result.json').is_file():
        raise ValueError('No completed model evaluation')
    shutil.copytree(source, destination, dirs_exist_ok=True)
    protocol = json.loads((root / 'detected_panel_protocol_v2.json').read_text())
    expected_cohorts = {f"detected_panel_{item['index']:03d}" for item in protocol['cohorts']
        if item['stage'] == 'development'}
    for variant in ['context', 'pair_only', 'baseline']:
        evaluation = destination / 'model_evaluation'
        evaluation_record = json.loads((evaluation / 'result.json').read_text())
        scorer_path = evaluation / f'{variant}_development_scorer_input.json'
        scorer_input = json.loads(scorer_path.read_text())
        cohorts = {cohort['cohort_id']: cohort for cohort in scorer_input['cohorts']}
        if set(cohorts) != expected_cohorts or len(cohorts) != len(scorer_input['cohorts']):
            raise ValueError('Scorer cohort roster differs from the frozen development split')
        files = {scorer_path}
        remote_bindings = {}
        if variant != 'baseline':
            predictions_path = evaluation / f'{variant}_development_predictions.json'
            raw_predictions = json.loads(predictions_path.read_text())
            files.add(predictions_path)
        for item in protocol['cohorts']:
            if item['stage'] != 'development':
                continue
            cohort_id = f"detected_panel_{item['index']:03d}"
            directory = destination / f"cohort_{item['index']:03d}"
            cohort_result = directory / 'result.json'
            record = json.loads(cohort_result.read_text())
            if record['item'] != item or record['state'] != 'COMPLETED':
                raise ValueError('Cohort result differs from the frozen protocol')
            files.add(cohort_result)
            for name in ['truth.json', 'identifiers.json', 'baseline_predictions.json']:
                local_path = directory / name
                remote_path = next(path for path in record['output_sha256']
                    if path.endswith(f"/cohort_{item['index']:03d}/{name}"))
                digest = hashlib.sha256(local_path.read_bytes()).hexdigest()
                if digest != record['output_sha256'][remote_path] or digest != evaluation_record['input_sha256'][remote_path]:
                    raise ValueError('Raw cohort table differs from simulation or evaluation bindings')
                remote_bindings[str(local_path.relative_to(root))] = {'remote_path': remote_path, 'sha256': digest}
                files.add(local_path)
            truth_path, ids_path = directory / 'truth.json', directory / 'identifiers.json'
            truth = json.loads(truth_path.read_text())
            if cohorts[cohort_id]['truth'] != truth:
                raise ValueError('Scorer truth differs from retained raw table')
            if variant == 'baseline':
                predictions_path = directory / 'baseline_predictions.json'
                expected = json.loads(predictions_path.read_text())
                files.add(predictions_path)
            else:
                expected = [{key: value for key, value in row.items() if key != 'cohort_id'}
                    for row in raw_predictions if row['cohort_id'] == cohort_id]
            if cohorts[cohort_id]['predictions'] != expected:
                raise ValueError('Scorer predictions differ from retained raw table')
            files.update([truth_path, ids_path])
        files.update([evaluation / 'model_calibration_freeze.json', evaluation / 'result.json',
            root / 'detected_panel_protocol_v2.json', root / 'detected_panel_source_manifest.json'])
        manifest = {'stage': 'DEVELOPMENT', 'variant': variant,
            'scorer_input': str(scorer_path.relative_to(root)),
            'remote_bound_files': remote_bindings,
            'files_sha256': {str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
                for path in sorted(files)}}
        output = root / f'detected_panel_scorer_inputs/{variant}.manifest.json'
        output.parent.mkdir(exist_ok=True)
        output.write_text(json.dumps(manifest, indent=2) + '\n')
        print(variant, scorer_path)


if __name__ == '__main__':
    main()
