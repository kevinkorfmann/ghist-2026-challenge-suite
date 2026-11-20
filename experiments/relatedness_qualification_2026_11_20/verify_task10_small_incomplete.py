"""Independently audit all eight cohorts, preserving failed inference as missing."""
import hashlib
import itertools
import json
import subprocess
from datetime import UTC, datetime
from pathlib import Path

from scripts.score_relatedness_cohorts import evaluate_cohorts
from scripts.score_relatedness_published_metric import evaluate_published_metric


base = Path('/tmp/task10-handoff-review-20261120')
raw_root = base / 'structural_small_export_v3_incomplete'
comparator = base / 'structural_small_comparator_v1'
load = lambda path: json.loads(path.read_text())
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
pair = lambda row: tuple(sorted((row['indID1'], row['indID2'])))
manifest = load(raw_root / 'manifest.json')
assert sha(raw_root / 'manifest.json') == 'c827bad87662be183460e2ce28bc51fc4cc7de2151a14f75d33d48c01ed9150f'
bindings = {}
for name, entry in manifest.items():
    assert sha(raw_root / name) == entry['sha256'], name
    bindings[entry['remote_path']] = entry['sha256']
bindings.update(load(comparator / 'result.json')['input_sha256'])
remote_root = '/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task10-development-20261120/'
for path in comparator.iterdir():
    bindings[remote_root + 'structural_small_comparator_v1/' + path.name] = sha(path)
failure = (raw_root / 'structural-small-10858090_0.out').read_text()
assert 'Initial NJ tree has no accepted root' in failure
reports = {}
for role, filename in [('structural', 'structural_predictions.json'), ('pair_only', 'pair_only_predictions.json'), ('baseline', 'baseline_predictions.json')]:
    data = load(comparator / f'{role}_scorer_input.json')
    assert len(data['cohorts']) == 8
    correct = total = missing = 0
    for index, cohort in enumerate(data['cohorts']):
        assert cohort['cohort_id'] == f'nested_{index:03d}'
        folder = raw_root / 'structural_qualification_v3' / f'cohort_{index:03d}'
        truth, ids = load(folder / 'truth.json'), load(folder / 'identifiers.json')
        expected = set(itertools.combinations(sorted(ids.values()), 2))
        assert len(ids) == len(set(ids.values()))
        assert truth == cohort['truth'] and len(truth) == len(expected)
        assert {pair(row) for row in truth} == expected
        if index == 0 and role == 'structural':
            predictions = []
        elif index == 0 and role == 'baseline':
            predictions = load(comparator / filename)
        else:
            predictions = load(folder / filename)
        assert predictions == cohort['predictions']
        assert len(predictions) == len({pair(row) for row in predictions})
        assert {pair(row) for row in predictions} <= expected
        if predictions:
            assert {pair(row) for row in predictions} == expected
        truth_labels = {pair(row): row['relation'] for row in truth}
        correct += sum(row['relation'] == truth_labels[pair(row)] for row in predictions)
        total += len(truth)
        missing += len(truth) - len(predictions)
    assert total == 1552
    expected_correct = {'structural': 1432, 'pair_only': 965, 'baseline': 730}[role]
    assert correct == expected_correct
    assert missing == (120 if role == 'structural' else 0)
    reports[role] = {'correct': correct, 'total': total, 'missing': missing,
        'scores': evaluate_cohorts(data['cohorts'], data['labels'], bootstrap_seed=20301120),
        'published': evaluate_published_metric(data['cohorts'], data['labels'])}
crop = load(base / 'crop_parity_result.json')
assert crop['state'] == 'PASS' and len(crop['cohorts']) == 8
for cohort in crop['cohorts']:
    assert cohort['crop']['records_excluded'] == 0
    assert len(set(cohort['comparison'].values())) == 1
bindings.update(crop['input_sha256'])
script = 'import pathlib,hashlib,json\nb=' + repr(bindings) + '\nfor p,s in b.items():\n assert hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()==s,p\nprint(json.dumps({"checked":len(b)}))\n'
remote = subprocess.run(['ssh', 'betty', 'python3 -'], input=script, text=True, capture_output=True, check=True)
report = {'utc': datetime.now(UTC).isoformat(), 'evidence_integrity': 'PASS',
    'strict_qualification': False, 'practical_scaling_gate': False,
    'reason': 'One cohort has 120 missing structural predictions; never exclude this cohort or claim complete coverage.',
    'raw_manifest_files': len(manifest), 'remote_hashes': json.loads(remote.stdout),
    'crop_identity_cohorts': 8, 'metrics': reports,
    'conditional_fallback': 'Not applied here; separate later policy requires prospective assessment.'}
Path(__file__).with_name('task10_small_incomplete_audit.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
print(json.dumps({k: v for k, v in report.items() if k != 'metrics'}))
for role, value in reports.items():
    print(role, value['correct'], value['total'], value['missing'])
