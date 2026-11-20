"""Independent Task11 output-only audit and deterministic four-column packaging."""
import ast
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
from itertools import combinations
import json
import math
from pathlib import Path
import sys

ROOT = Path('/private/tmp/ghist-final-submissions-20261120/task11')
INPUTS = ROOT / 'audit_inputs'
OUTPUTS = INPUTS / 'outputs'
SOURCE = Path('/tmp/ghist-task10-final-independent-20261120/source')
CANDIDATE_SHA = 'f2c1391a5866567a3bc10423bea25ecd537b86f9d3cf819fcd72798bb2abd833'
ORIGINAL_SHA = '5bbdd1a850f65308ed2376a159cbe63c7540e76eab9e91bb77b75c127a2fc553'
MAPPED_SHA = '24df5742918baa0c3e7a5725e87dc904a31429f66b320be455f26115003618cd'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def write_new(path, record):
    with path.open('x') as stream:
        stream.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + '\n')


def main():
    candidate_path = INPUTS / 'pragmatic_corrected_context_candidate_v1.json'
    assert sha(candidate_path) == CANDIDATE_SHA
    candidate, adapter = read(candidate_path), read(OUTPUTS / 'adapter_audit.json')
    assert candidate['selected_role'] == adapter['selected_role'] == 'context'
    assert adapter['candidate_sha256'] == CANDIDATE_SHA and adapter['completed'] is True
    assert adapter['job_id'] == '10859650' and adapter['truth_used'] is False
    assert adapter['input_vcf_sha256'] == MAPPED_SHA
    assert adapter['protocol_sha256'] == candidate['protocol_sha256']
    for role in ['base', 'context']:
        assert adapter[f'{role}_model_sha256'] == candidate['bound_inputs'][candidate[f'{role}_model']]
    assert adapter['adapter_sha256'] == sha(INPUTS / 'prepare_pragmatic_corrected_v1.py')
    for name in ['prepare_pragmatic_corrected_v1.py', 'pragmatic_corrected_bindings_v1.py']:
        remote = '/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task11-development-20261120/' + name
        assert sha(INPUTS / name) == candidate['bound_inputs'][remote]
    remote_bindings = read(INPUTS / 'remote_binding_hashes.json')
    local_bindings = {}
    for name, expected in remote_bindings.items():
        if '/source_v5/src/' in name:
            local = SOURCE / name.split('/source_v5/src/')[1]
            assert expected == candidate['source_sha256']['src/' + str(local.relative_to(SOURCE))]
        else:
            local = INPUTS / Path(name).name
        assert sha(local) == expected, name
        local_bindings[str(local)] = expected
    remote_outputs = read(OUTPUTS / 'remote_compact_hashes.json')
    assert remote_outputs['all_compact_outputs_match'] is True
    for name, expected in remote_outputs['compact_output_sha256'].items():
        assert sha(OUTPUTS / name) == expected, name
        if name != 'adapter_audit.json':
            assert adapter['output_sha256'][name] == expected
    assert set(adapter['output_sha256']) == (
        set(remote_outputs['compact_output_sha256']) - {'adapter_audit.json'}
        | set(remote_outputs['vcf_hashes_deferred_to_root']))

    preflight, provenance = read(INPUTS / 'vcf_preflight.json'), read(INPUTS / 'input_provenance.json')
    assert preflight['input_sha256'] == provenance['sha256'] == ORIGINAL_SHA
    original = preflight['sample_ids']
    assert preflight['sample_count'] == len(original) == len(set(original)) == 68
    roster = read(OUTPUTS / 'roster.json')
    assert roster['samples'] == original and len(set(roster['internal_samples'])) == 68
    assert dict(zip(roster['internal_samples'], original)) == roster['internal_to_original']
    family_record = read(OUTPUTS / 'families.json')
    families = family_record['families']
    assert len(families) == family_record['supplied_family_count'] == 8
    flattened = [person for family in families for person in family]
    assert len(flattened) == len(set(flattened)) == 68 and set(flattened) == set(original)
    clusters = read(OUTPUTS / 'clusters.json')['families']
    assert [[roster['internal_to_original'][person] for person in family] for family in clusters] == families
    expected = {tuple(sorted(pair)) for family in families for pair in combinations(family, 2)}
    assert len(expected) == 265 == adapter['pairs'] and adapter['samples'] == 68

    node = next(node for node in ast.parse((SOURCE / 'ghist_relatedness/labels.py').read_text()).body
                if isinstance(node, ast.Assign)
                and any(isinstance(target, ast.Name) and target.id == 'DOG_LABELS' for target in node.targets))
    labels = ast.literal_eval(node.value)
    meioses = {label: up + down for (up, down, _), label in labels.items()}
    probability_checks = {}
    rows = {}
    for filename in ['base_predictions.json', 'context_predictions.json', 'predictions.json']:
        records = read(OUTPUTS / filename)
        rows[filename] = records
        assert len(records) == len(expected)
        counts = Counter()
        max_sum_error = 0
        for row in records:
            pair = tuple(sorted((row['indID1'], row['indID2'])))
            counts[pair] += 1
            assert row['indID1'] != row['indID2']
            assert row['indID1'] in families[row['family']] and row['indID2'] in families[row['family']]
            probabilities = row['probabilities']
            assert set(probabilities) == set(meioses)
            assert all(not isinstance(p, bool) and isinstance(p, (int, float))
                       and math.isfinite(p) and 0 <= p <= 1 for p in probabilities.values())
            error = abs(sum(probabilities.values()) - 1)
            assert error < 1e-10
            max_sum_error = max(max_sum_error, error)
            assert row['relation'] in probabilities
            assert math.isfinite(row['support']) and math.isclose(row['support'], probabilities[row['relation']], rel_tol=1e-12, abs_tol=1e-12)
        assert set(counts) == expected and set(counts.values()) == {1}
        probability_checks[filename] = {'rows': len(records), 'all_finite_nonnegative_normalized': True,
                                      'maximum_probability_sum_error': max_sum_error,
                                      'support_matches_selected_label_probability': True}
    assert rows['predictions.json'] == rows['context_predictions.json']
    with (OUTPUTS / 'task11_candidate.csv').open(newline='') as stream:
        raw = list(csv.reader(stream, strict=True))
    assert raw.pop(0) == ['indID1', 'indID2', 'relation']
    triples = [[row['indID1'], row['indID2'], row['relation']] for row in rows['predictions.json']]
    assert raw == triples

    sys.path.insert(0, str(SOURCE))
    from ghist_submission.validate import validate_submission, within_family_pairs
    import ghist_submission.validate as validator
    assert Path(validator.__file__).resolve() == (SOURCE / 'ghist_submission/validate.py').resolve()
    assert within_family_pairs(families, original) == expected
    raw_validation = validate_submission(OUTPUTS / 'task11_candidate.csv', 11, 'final', expected)
    final = ROOT / 'task11_final_candidate.csv'
    with final.open('x', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['indID1', 'indID2', 'relation', 'meioses_count'])
        writer.writerows([*row, meioses[row[2]]] for row in raw)
    with final.open(newline='') as stream:
        packaged = list(csv.reader(stream, strict=True))
    assert packaged.pop(0) == ['indID1', 'indID2', 'relation', 'meioses_count']
    assert [row[:3] for row in packaged] == raw
    assert all(row[3] == str(meioses[row[2]]) for row in packaged)
    validation = validate_submission(final, 11, 'final', expected)
    assert validation == raw_validation == read(OUTPUTS / 'validation.json')

    report = {'state': 'PASS_INDEPENDENT_OUTPUT_AUDIT_AND_DETERMINISTIC_PACKAGING',
        'verified_at_utc': datetime.now(timezone.utc).isoformat(),
        'scope': 'Existing frozen outputs and original preflight roster only; no actual VCF, truth, tuning, inference, or upload.',
        'candidate_local_path': str(final), 'candidate_sha256': sha(final),
        'frozen_model_candidate_sha256': CANDIDATE_SHA,
        'original_input_sha256_from_root_preflight': ORIGINAL_SHA,
        'mapped_input_sha256_from_frozen_completed_adapter': MAPPED_SHA,
        'original_unique_samples': 68, 'families': 8, 'family_sizes': [len(f) for f in families],
        'within_family_pairs': 265, 'all_within_inferred_family_pairs_exactly_once': True,
        'context_predictions_equal_selected_predictions': True,
        'raw_and_packaged_first_three_columns_identical_to_context_labels_in_order': True,
        'meioses_mapping': meioses, 'meioses_rule': 'Frozen DOG_LABELS tuple up + down',
        'probability_checks': probability_checks, 'frozen_official_schema_validation': validation,
        'compact_output_hash_count': len(remote_outputs['compact_output_sha256']),
        'all_compact_outputs_match_adapter_and_current_remote_sha256': True,
        'vcf_hashes_deferred_to_root': remote_outputs['vcf_hashes_deferred_to_root'],
        'source_binding_local_remote_sha256': local_bindings,
        'compact_artifact_sha256': {name: sha(OUTPUTS / name) for name in remote_outputs['compact_output_sha256']},
        'packaging_audit_source_sha256': sha(Path(__file__)),
        'strict_qualification_passed': False, 'submission_performed': False,
        'limitation': 'Valid packaging and normalized probabilities do not establish inference accuracy; root retains source/model/full-output rehash and submission responsibility.'}
    provenance_record = {'scope': 'Deterministic meiosis column only; all predictions unchanged',
        'created_at_utc': report['verified_at_utc'], 'input_csv_sha256': sha(OUTPUTS / 'task11_candidate.csv'),
        'selected_predictions_sha256': sha(OUTPUTS / 'predictions.json'),
        'context_predictions_sha256': sha(OUTPUTS / 'context_predictions.json'),
        'preflight_sha256': sha(INPUTS / 'vcf_preflight.json'), 'adapter_audit_sha256': sha(OUTPUTS / 'adapter_audit.json'),
        'frozen_model_candidate_sha256': CANDIDATE_SHA,
        'frozen_labels_sha256': sha(SOURCE / 'ghist_relatedness/labels.py'),
        'frozen_validator_sha256': sha(SOURCE / 'ghist_submission/validate.py'),
        'output_path': str(final), 'output_sha256': sha(final), 'meioses_mapping': meioses,
        'validation': validation, 'packaging_source_sha256': sha(Path(__file__)), 'submission_performed': False}
    destination = Path(__file__).parent
    write_new(destination / 'task11_final_packaging_provenance.json', provenance_record)
    write_new(destination / 'task11_final_packaging_audit.json', report)
    print(json.dumps({'state': report['state'], 'candidate_local_path': str(final),
        'candidate_sha256': sha(final), 'audit_sha256': sha(destination / 'task11_final_packaging_audit.json')}))


if __name__ == '__main__':
    main()
