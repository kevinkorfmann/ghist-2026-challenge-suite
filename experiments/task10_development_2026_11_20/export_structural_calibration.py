"""Bind calibration scorer rows to original truth, raw inference and global-q rule."""
import hashlib
import itertools
import json
import shutil
from pathlib import Path


def main():
    root = Path(__file__).resolve().parent
    source = Path('/tmp/task10-evidence-20261120/structural_calibration_v3')
    cache = Path('/tmp/task10-evidence-20261120/detected_panel_v2')
    destination = root / 'structural_calibration_evidence'
    if destination.exists():
        raise ValueError('Calibration export already exists')
    shutil.copytree(source, destination)
    def digest(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()
    confidence = json.loads((destination / 'confidence_freeze.json').read_text())
    scorer = json.loads((destination / 'scorer_input.json').read_text())
    files = set(destination.rglob('*.json'))
    correct, pairs = 0, 0
    for cohort in scorer['cohorts']:
        index = int(cohort['cohort_id'].rsplit('_', 1)[1])
        paths = {}
        for name in ['truth.json', 'identifiers.json']:
            directory = destination / f'cohort_{index:03d}'
            directory.mkdir(exist_ok=True)
            path = paths[name] = directory / name
            shutil.copyfile(cache / f'cohort_{index:03d}' / name, path)
            files.add(path)
        remote_truth = next(path for path in confidence['input_sha256'] if path.endswith(f'/cohort_{index:03d}/truth.json'))
        if digest(paths['truth.json']) != confidence['input_sha256'][remote_truth]:
            raise ValueError('Original calibration truth binding differs')
        if cohort['truth'] != json.loads(paths['truth.json'].read_text()):
            raise ValueError('Scorer truth differs from original raw table')
        samples = sorted(json.loads(paths['identifiers.json'].read_text()).values())
        expected_pairs = set(itertools.combinations(samples, 2))
        observed = {tuple(sorted([row['indID1'], row['indID2']])) for row in cohort['truth']}
        if observed != expected_pairs or len(cohort['truth']) != len(expected_pairs):
            raise ValueError('Incomplete calibration roster pair coverage')
        raw = json.loads((destination / f'cohort_{index:03d}_inference.json').read_text())['predictions']
        q = confidence['confidence']
        expected = [{**row, 'probabilities': {label: q if row['relation'] == label else (1 - q) / 13
            for label in scorer['labels']}} for row in raw]
        if expected != cohort['predictions']:
            raise ValueError('Scorer differs from raw inference and frozen q rule')
        keyed = {(row['indID1'], row['indID2']): row['relation'] for row in cohort['truth']}
        correct += sum(keyed[(row['indID1'], row['indID2'])] == row['relation'] for row in raw)
        pairs += len(raw)
    if pairs != 623 or correct != confidence['correct'] or (correct + 1) / (pairs + 2) != confidence['confidence']:
        raise ValueError('Global correctness rule differs from complete calibration rows')
    for name in ['structural_qualification_protocol_v3.json', 'structural_confidence.py',
            'calibrate_structural_candidate.py', 'structural_pair_candidate.py',
            'integer_edge_decoder.py', 'forest_distance_decoder.py']:
        files.add(root / name)
    output = root / 'structural_calibration_manifest.json'
    output.write_text(json.dumps({'stage': 'CALIBRATION_NOT_INDEPENDENT_QUALIFICATION',
        'scorer_input': 'structural_calibration_evidence/scorer_input.json',
        'confidence_interpretation': 'One global correctness estimate; zero individual resolution; no root posterior',
        'files_sha256': {str(path.relative_to(root)): digest(path) for path in sorted(files)}}, indent=2) + '\n')


if __name__ == '__main__':
    main()
