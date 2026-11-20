"""Independent binding, raw-pair and constant-confidence calibration audit."""
import hashlib
import itertools
import json
from datetime import UTC, datetime
from pathlib import Path

from scripts.score_relatedness_cohorts import evaluate_cohorts
from scripts.score_relatedness_published_metric import evaluate_published_metric


base = Path('/tmp/task10-review-v2-20261120')
manifest_path = base / 'structural_calibration_manifest.json'
manifest_bytes = manifest_path.read_bytes()
manifest = json.loads(manifest_bytes)
assert manifest['stage'] == 'CALIBRATION_NOT_INDEPENDENT_QUALIFICATION'
for name, expected in manifest['files_sha256'].items():
    assert not Path(name).is_absolute() and '..' not in Path(name).parts
    assert hashlib.sha256((base / name).read_bytes()).hexdigest() == expected, name
root = base / 'structural_calibration_evidence'
freeze = json.loads((root / 'confidence_freeze.json').read_text())
raw = json.loads((base / manifest['scorer_input']).read_text())
assert len(raw['labels']) == len(set(raw['labels'])) == 14
assert {x['cohort_id'] for x in raw['cohorts']} == {f'detected_panel_{i:03d}' for i in range(32, 40)}
correct = count = 0
for cohort in raw['cohorts']:
    i = int(cohort['cohort_id'].rsplit('_', 1)[1])
    truth = json.loads((root / f'cohort_{i:03d}/truth.json').read_text())
    ids = json.loads((root / f'cohort_{i:03d}/identifiers.json').read_text())
    inferred = json.loads((root / f'cohort_{i:03d}_inference.json').read_text())['predictions']
    assert truth == cohort['truth']
    assert len(ids) == len(set(ids.values()))
    expected_pairs = set(itertools.combinations(sorted(ids.values()), 2))
    pair = lambda row: tuple(sorted((row['indID1'], row['indID2'])))
    for rows in (truth, inferred, cohort['predictions']):
        assert len(rows) == len(expected_pairs)
        assert {pair(row) for row in rows} == expected_pairs
    assert inferred == [{k: v for k, v in row.items() if k != 'probabilities'} for row in cohort['predictions']]
    lookup = {pair(row): row['relation'] for row in truth}
    for row in cohort['predictions']:
        assert set(row['probabilities']) == set(raw['labels'])
        for label, value in row['probabilities'].items():
            expected = freeze['confidence'] if label == row['relation'] else (1 - freeze['confidence']) / 13
            assert abs(value - expected) < 1e-14
        correct += row['relation'] == lookup[pair(row)]
        count += 1
assert correct == freeze['correct'] == 483
assert count == freeze['pairs'] == 623
assert freeze['confidence'] == (correct + 1) / (count + 2) == .7744
assert freeze['individual_resolution'] == 0 and freeze['root_posterior'] is False
scores = evaluate_cohorts(raw['cohorts'], raw['labels'], bootstrap_seed=20301120)
published = evaluate_published_metric(raw['cohorts'], raw['labels'])
assert scores['pooled']['correct'] == correct
for name, expected in manifest['files_sha256'].items():
    assert hashlib.sha256((base / name).read_bytes()).hexdigest() == expected, name
assert manifest_path.read_bytes() == manifest_bytes
report = {
    'verified_at_utc': datetime.now(UTC).isoformat(),
    'scope': 'Calibration only, not independent qualification; all original roster, truth and prediction rows verified.',
    'manifest_sha256': hashlib.sha256(manifest_bytes).hexdigest(),
    'bound_file_count': len(manifest['files_sha256']),
    'scores': scores, 'published_metric': published,
    'correct': correct, 'pairs': count, 'q': freeze['confidence'],
    'individual_resolution': 0, 'root_posterior': False,
    'remote_bindings_to_verify': freeze['input_sha256'],
}
out = Path(__file__).with_name('task10_structural_calibration_audit.json')
out.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
print(json.dumps({k: report[k] for k in ['correct', 'pairs', 'q', 'bound_file_count']}))
