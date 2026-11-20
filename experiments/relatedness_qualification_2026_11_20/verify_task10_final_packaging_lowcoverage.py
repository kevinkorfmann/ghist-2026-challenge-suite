"""Independent output-only Task10 packaging audit; no VCF/model/truth access."""
import ast
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
from itertools import combinations
import json
from pathlib import Path
import sys


LOCAL = Path('/tmp/ghist-task10-final-independent-20261120')
REMOTE = '/vast/projects/smathi/cohort/korfmann/ghist2026/final_runs/task10-root-20261120-2028'
SOURCE = '/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task10-development-20261120/detected-panel-source'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def main():
    packaging = read(LOCAL / 'task10_final_packaging_provenance.json')
    remote_hashes = read(LOCAL / 'remote_current_hashes.json')
    assert remote_hashes['all_current_remote_packaging_bindings_match'] is True
    mapping = {
        f'{REMOTE}/task10_final_candidate.csv': LOCAL / 'task10_final_candidate.csv',
        f'{REMOTE}/task10_final_packaging_provenance.json': LOCAL / 'task10_final_packaging_provenance.json',
        f'{REMOTE}/vcf_preflight.json': LOCAL / 'vcf_preflight.json',
        f'{SOURCE}/source_manifest.json': LOCAL / 'source/source_manifest.json',
    }
    for name in ['candidate.csv', 'structural_predictions.json', 'result.json', 'initial_provenance.json']:
        mapping[f'{REMOTE}/pragmatic_fallback_v1_run/{name}'] = LOCAL / 'raw' / name
    for name in ['ghist_submission/__init__.py', 'ghist_submission/validate.py',
                 'ghist_relatedness/__init__.py', 'ghist_relatedness/labels.py']:
        mapping[f'{SOURCE}/src/{name}'] = LOCAL / 'source' / name
    for remote, local in mapping.items():
        assert sha(local) == remote_hashes['sha256'][remote], remote
        if remote in packaging['input_sha256']:
            assert sha(local) == packaging['input_sha256'][remote], remote
    manifest = read(LOCAL / 'source/source_manifest.json')
    for remote, local in mapping.items():
        if remote.startswith(SOURCE + '/src/'):
            assert sha(local) == manifest[remote], remote

    preflight = read(LOCAL / 'vcf_preflight.json')
    roster = preflight['sample_ids']
    assert preflight['sample_count'] == len(roster) == len(set(roster)) == 263
    expected_pairs = set(combinations(sorted(roster), 2))
    assert len(expected_pairs) == 34453
    initial = read(LOCAL / 'raw/initial_provenance.json')
    result = read(LOCAL / 'raw/result.json')
    assert initial['input_sha256'] == preflight['input_sha256']
    assert result['input_and_model_confidence_sha256'][preflight['input_path']] == preflight['input_sha256']
    assert initial['protocol_sha256'] == result['protocol_sha256']
    assert sha(LOCAL / 'raw/candidate.csv') == result['candidate_sha256']
    for name in ['candidate.csv', 'structural_predictions.json', 'initial_provenance.json']:
        assert sha(LOCAL / 'raw' / name) == result['output_sha256'][f'{REMOTE}/pragmatic_fallback_v1_run/{name}']

    label_source = LOCAL / 'source/ghist_relatedness/labels.py'
    label_assignments = [node for node in ast.parse(label_source.read_text()).body
                         if isinstance(node, ast.Assign)
                         and any(isinstance(target, ast.Name) and target.id == 'HUMAN_LABELS' for target in node.targets)]
    assert len(label_assignments) == 1
    human_labels = ast.literal_eval(label_assignments[0].value)
    meioses = {}
    for (up, down, _), label in human_labels.items():
        assert label not in meioses or meioses[label] == up + down
        meioses[label] = up + down
    assert meioses == packaging['label_to_meioses_count']

    with (LOCAL / 'raw/candidate.csv').open(newline='') as stream:
        raw = list(csv.reader(stream, strict=True))
    with (LOCAL / 'task10_final_candidate.csv').open(newline='') as stream:
        packaged = list(csv.reader(stream, strict=True))
    assert raw.pop(0) == ['indID1', 'indID2', 'relation']
    assert packaged.pop(0) == ['indID1', 'indID2', 'relation', 'meioses_count']
    predictions = read(LOCAL / 'raw/structural_predictions.json')
    assert len(raw) == len(packaged) == len(predictions) == 34453
    counts = Counter()
    for original, final, prediction in zip(raw, packaged, predictions):
        assert len(original) == 3 and len(final) == 4
        assert final[:3] == original == [prediction['indID1'], prediction['indID2'], prediction['relation']]
        left, right, label = original
        assert left != right and left in roster and right in roster
        assert final[3] == str(meioses[label])
        counts[tuple(sorted((left, right)))] += 1
    assert set(counts) == expected_pairs and set(counts.values()) == {1}

    sys.path.insert(0, str(LOCAL / 'source'))
    from ghist_submission.validate import validate_relatedness
    import ghist_submission.validate as validator
    assert Path(validator.__file__).resolve() == (LOCAL / 'source/ghist_submission/validate.py').resolve()
    validation = validate_relatedness(LOCAL / 'task10_final_candidate.csv', 10, expected_pairs=expected_pairs)
    assert validation == {'pairs': 34453, 'individuals': 263, 'input_pair_coverage_verified': True}
    candidate_sha = sha(LOCAL / 'task10_final_candidate.csv')
    assert candidate_sha == packaging['output_sha256']
    report = {
        'state': 'PASS_INDEPENDENT_OUTPUT_ONLY_PACKAGING_AUDIT',
        'verified_at_utc': datetime.now(timezone.utc).isoformat(),
        'auditor': 'lowcoverage_readiness, independent of Task10 inference and packaging',
        'scope': 'Compact frozen outputs and original preflight roster only; no VCF, truth, model fitting, tuning, or upload.',
        'candidate_local_path': str(LOCAL / 'task10_final_candidate.csv'),
        'candidate_remote_path': packaging['output_path'], 'candidate_sha256': candidate_sha,
        'original_unique_roster_count': 263, 'unordered_pair_count': 34453,
        'every_original_roster_pair_exactly_once': True,
        'raw_prediction_to_original_and_packaged_first_three_columns_exact_in_order': True,
        'meioses_equal_independently_parsed_frozen_human_labels_up_plus_down': True,
        'frozen_schema_validator_result': validation,
        'frozen_schema_validator_sha256': sha(LOCAL / 'source/ghist_submission/validate.py'),
        'frozen_human_labels_sha256': sha(label_source),
        'frozen_source_manifest_sha256': sha(LOCAL / 'source/source_manifest.json'),
        'original_input_hash_consistent_across_preflight_initial_and_result': True,
        'local_remote_sha256_match_count': len(mapping),
        'all_current_remote_packaging_bindings_match': True,
        'bound_local_artifacts': {str(path): sha(path) for path in mapping.values()},
        'remote_hash_audit_sha256': sha(LOCAL / 'remote_current_hashes.json'),
        'audit_source_sha256': sha(Path(__file__)),
        'scientific_qualification': 'Unchanged; packaging verification does not validate prediction accuracy.',
        'submission_performed': False,
    }
    output = Path(__file__).with_name('task10_final_packaging_audit.json')
    with output.open('x') as stream:
        stream.write(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'state': report['state'], 'candidate_local_path': report['candidate_local_path'],
                      'candidate_sha256': candidate_sha, 'audit_path': str(output), 'audit_sha256': sha(output)}))


if __name__ == '__main__':
    main()
