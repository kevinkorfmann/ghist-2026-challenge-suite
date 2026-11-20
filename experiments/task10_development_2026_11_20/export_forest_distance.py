"""Export supplementary decoder rows with raw tables and original truth bindings."""
import hashlib
import json
import shutil
from pathlib import Path


def main():
    root = Path(__file__).resolve().parent
    source = Path('/tmp/task10-evidence-20261120/forest_distance_decoder_v1')
    destination = root / 'forest_distance_evidence'
    shutil.copytree(source, destination, dirs_exist_ok=True)
    result = json.loads((destination / 'result.json').read_text())
    if result.get('state') == 'SKIPPED_PRIMARY_V2_PASSED':
        print('Skipped; primary v2 passed')
        return
    for variant in ['context', 'pair_only']:
        scorer_path = destination / f'{variant}_scorer_input.json'
        scorer = json.loads(scorer_path.read_text())
        raw_path = destination / f'{variant}_predictions.json'
        raw = json.loads(raw_path.read_text())
        primary_manifest_path = root / f'detected_panel_scorer_inputs/{variant}.manifest.json'
        primary = json.loads(primary_manifest_path.read_text())
        files = {scorer_path, raw_path, destination / f'{variant}_trees.json', destination / 'result.json',
            primary_manifest_path, root / 'forest_distance_protocol.json', root / 'forest_distance_decoder.py'}
        for name, expected in primary['files_sha256'].items():
            path = root / name
            if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
                raise ValueError('Primary retained evidence changed')
        for cohort in scorer['cohorts']:
            index = int(cohort['cohort_id'].rsplit('_', 1)[1])
            directory = root / f'detected_panel_evidence/cohort_{index:03d}'
            truth_path, ids_path = directory / 'truth.json', directory / 'identifiers.json'
            if cohort['truth'] != json.loads(truth_path.read_text()):
                raise ValueError('Decoder truth differs from retained raw truth')
            expected = [{key: value for key, value in row.items() if key != 'cohort_id'}
                for row in raw if row['cohort_id'] == cohort['cohort_id']]
            if expected != cohort['predictions']:
                raise ValueError('Decoder scorer rows differ from retained raw predictions')
            files.update([truth_path, ids_path])
        output = root / f'forest_distance_scorer_inputs/{variant}.manifest.json'
        output.parent.mkdir(exist_ok=True)
        output.write_text(json.dumps({'stage': 'REUSED_DEVELOPMENT_DIAGNOSTIC', 'variant': variant,
            'scorer_input': str(scorer_path.relative_to(root)), 'uncertainty_calibrated': False,
            'files_sha256': {str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
                for path in sorted(files)}}, indent=2) + '\n')
        print(variant, scorer_path)


if __name__ == '__main__':
    main()
