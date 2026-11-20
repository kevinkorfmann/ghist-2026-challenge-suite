"""Independently bind and score the retained Task 10 development comparisons."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from itertools import combinations
from pathlib import Path


ROOT = Path('/Users/kevinkorfmann/Documents/ChatGPT/ghist')
BASE = ROOT / 'experiments/task10_development_2026_11_20'
OUT = ROOT / 'experiments/relatedness_qualification_2026_11_20/task10_development'
SCORER = Path('/private/tmp/ghist-root-qualification-20261120/score_relatedness_cohorts.py')
LABEL_SOURCE = Path('/private/tmp/ghist-task10-development-20261120/src/ghist_relatedness/labels.py')


def load_module(name, path):
    specification = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical(rows):
    return sorted(rows, key=lambda row: tuple(sorted((row['indID1'], row['indID2']))))


def main():
    scorer = load_module('independent_scorer', SCORER)
    labels = load_module('verified_labels', LABEL_SOURCE)
    OUT.mkdir(parents=True, exist_ok=False)
    summary = {}
    for name in ('raw_haplotype', 'diploid_state', 'exact_diploid_oracle', 'whole_family_bonsai', 'NJ_5', 'NJ_10'):
        manifest_path = BASE / f'scorer_inputs/{name}.manifest.json'
        manifest_bytes = manifest_path.read_bytes()
        manifest = json.loads(manifest_bytes)
        assert manifest['stage'] == 'DEVELOPMENT'
        retained = {}
        for relative, expected in manifest['files_sha256'].items():
            path = (BASE / relative).resolve()
            assert path.is_relative_to(BASE.resolve())
            data = path.read_bytes()
            assert hashlib.sha256(data).hexdigest() == expected, relative
            retained[relative] = json.loads(data)
        data = retained[f'scorer_inputs/{name}.json']
        assert set(data['labels']) == set(labels.HUMAN_LABELS.values())
        source_predictions = {p: r for p, r in retained.items()
                              if not p.startswith('scorer_inputs/') and not p.endswith('/truth.json')}
        for index, cohort in enumerate(data['cohorts']):
            truth_path = f'evidence/development_cohort_{index}/truth.json'
            assert canonical(cohort['truth']) == canonical(retained[truth_path])
            matches = [p for p, r in source_predictions.items()
                       if canonical(cohort['predictions']) == canonical(r)]
            assert len(matches) == 1, (name, index, matches)
            roster = {r[k] for r in cohort['truth'] for k in ('indID1', 'indID2')}
            observed = {tuple(sorted((r['indID1'], r['indID2']))) for r in cohort['truth']}
            assert observed == set(combinations(sorted(roster), 2))
        report = scorer.evaluate_cohorts(data['cohorts'], data['labels'])
        report.update(stage='DEVELOPMENT', source_manifest_sha256=hashlib.sha256(manifest_bytes).hexdigest(),
                      source_hashes=manifest['files_sha256'], scorer_sha256=digest(SCORER),
                      label_source_sha256=digest(LABEL_SOURCE), raw_rows_verified=True,
                      complete_truth_rosters_verified=True)
        for relative, expected in manifest['files_sha256'].items():
            assert digest(BASE / relative) == expected, relative
        assert manifest_path.read_bytes() == manifest_bytes
        (OUT / f'{name}.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
        summary[name] = {'pairs': report['pooled']['pairs'], 'correct': report['pooled']['correct'],
                         'accuracy': report['pooled']['exact_accuracy'],
                         'balanced_accuracy': report['pooled']['balanced_accuracy'],
                         'minimum_label_recall': min(r['recall'] for r in report['pooled']['per_label'].values()),
                         'coverage_passed': report['pair_coverage_passed'],
                         'all_labels_represented': report['pooled']['all_labels_represented'],
                         'ece': report['calibration']['ece'] if report['calibration'] else None,
                         'bootstrap_status': report['bootstrap']['status']}
    (OUT / 'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
