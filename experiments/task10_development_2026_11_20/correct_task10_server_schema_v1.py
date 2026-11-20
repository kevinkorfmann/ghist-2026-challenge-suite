"""Apply only the root-reported server allowlist; preserve the invalid submission."""

from collections import Counter
from pathlib import Path
import csv
import datetime
import hashlib
import itertools
import json
import os
import tempfile


ALLOWED = (
    'full_1st_cousin', 'full_2nd_cousin', 'siblings', 'avuncular',
    'full_1st_cousin_once_removed', 'full_2nd_cousin_once_removed',
    'full_3rd_cousin', 'full_2nd_cousin_twice_removed',
    'full_3rd_cousin_once_removed', 'full_1st_cousin_twice_removed',
    'full_4th_cousin', 'avuncular_once_removed',
)
MAPPING = {
    'full_1st_cousin_thrice_removed': 'full_2nd_cousin_once_removed',
    'avuncular_twice_removed': 'full_1st_cousin_once_removed',
}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parent = Path('/vast/projects/smathi/cohort/korfmann/ghist2026/final_runs/task10-root-20261120-2028')
    original = parent / 'task10_final_candidate.csv'
    packaging = parent / 'task10_final_packaging_provenance.json'
    preflight = parent / 'vcf_preflight.json'
    output = parent / 'task10_final_schema_corrected.csv'
    report = parent / 'task10_final_schema_correction_provenance.json'
    if output.exists() or report.exists():
        raise FileExistsError('Schema correction outputs already exist')
    bound = {str(p): sha256(p) for p in [original, packaging, preflight, Path(__file__).resolve()]}
    if bound[str(original)] != 'ee9bace7bfc086281be0b243c5b440db9542277466b4722a14cb32523ef059be':
        raise ValueError('Unexpected original submission bytes')
    if bound[str(packaging)] != '1643526d98cb9d270f53f493c66271500cff94164a718049123f99a6e16bb931':
        raise ValueError('Unexpected original packaging provenance')
    counts = json.loads(packaging.read_text())['label_to_meioses_count']
    if len(ALLOWED) != 12 or len(set(ALLOWED)) != 12 or set(counts) - set(ALLOWED) != set(MAPPING):
        raise ValueError('Unexpected allowlist or unsupported labels')
    for old, new in MAPPING.items():
        equivalents = [label for label in ALLOWED if counts[label] == counts[old]]
        if equivalents != [new]:
            raise ValueError('Replacement is not the unique allowed same-meiosis label')
    roster = json.loads(preflight.read_text())['sample_ids']
    if len(roster) != len(set(roster)) or len(roster) != 263:
        raise ValueError('Expected the full original 263-person roster')
    expected = set(itertools.combinations(sorted(roster), 2))
    with original.open(newline='') as stream:
        reader = csv.DictReader(stream, strict=True)
        fields = reader.fieldnames
        if fields != ['indID1', 'indID2', 'relation', 'meioses_count']:
            raise ValueError('Unexpected original columns')
        rows = list(reader)
    pairs = [tuple(sorted([row['indID1'], row['indID2']])) for row in rows]
    if len(rows) != 34453 or len(set(pairs)) != len(rows) or set(pairs) != expected:
        raise ValueError('Original pair coverage mismatch')
    replacements = Counter()
    corrected = []
    for row in rows:
        if set(row) != set(fields) or int(row['meioses_count']) != counts[row['relation']]:
            raise ValueError('Original label/meiosis mismatch')
        new = MAPPING.get(row['relation'], row['relation'])
        if new != row['relation']:
            replacements[row['relation']] += 1
        if new not in ALLOWED or counts[new] != int(row['meioses_count']):
            raise ValueError('Corrected row violates server allowlist or meiosis count')
        corrected.append(row | {'relation': new})
    descriptor, name = tempfile.mkstemp(prefix='.task10-schema-', suffix='.csv', dir=parent)
    temporary = Path(name)
    try:
        with os.fdopen(descriptor, 'w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows(corrected)
        with temporary.open(newline='') as stream:
            check = list(csv.DictReader(stream, strict=True))
        if check != corrected:
            raise ValueError('CSV round-trip mismatch')
        for before, after in zip(rows, check):
            if any(before[key] != after[key] for key in ['indID1', 'indID2', 'meioses_count']):
                raise ValueError('Pair, order or meiosis count changed')
            if after['relation'] != MAPPING.get(before['relation'], before['relation']):
                raise ValueError('Unapproved relationship change')
        if any(sha256(p) != d for p, d in bound.items()):
            raise ValueError('Original artifact changed during correction')
        os.link(temporary, output)
    finally:
        temporary.unlink(missing_ok=True)
    provenance = {
        'state': 'SERVER_ALLOWLIST_SCHEMA_CORRECTED_PENDING_ROOT_REVIEW',
        'created_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'server_evidence': {'source': 'Root-reported rejection email', 'email_id': '1a1fad9451a8e239',
                            'time_utc': '2026-11-20T21:25', 'submission_id': '9785356', 'queue_id': '9619674',
                            'invalid_submission_not_scored': True},
        'server_allowed_labels': list(ALLOWED),
        'mapping': MAPPING,
        'mapping_rule': 'Unique server-accepted label with identical frozen meiosis count; this merges different pedigree labels, not biological equivalence',
        'replacement_counts': dict(replacements),
        'pairs': len(rows), 'people': len(roster),
        'all_pair_ids_order_and_meioses_preserved': True,
        'all_other_relations_preserved': True,
        'no_model_or_threshold_changes': True,
        'no_truth_or_scientific_score_feedback_used': True,
        'input_sha256_before_after': bound,
        'output_path': str(output), 'output_sha256': sha256(output),
        'scientific_qualification': False, 'submission_performed': False,
    }
    with report.open('x') as stream:
        json.dump(provenance, stream, indent=2)
        stream.write('\n')
    print(json.dumps({'output': str(output), 'sha256': sha256(output), 'provenance_sha256': sha256(report),
                      'replacement_counts': dict(replacements), 'pairs': len(rows)}))


if __name__ == '__main__':
    main()
