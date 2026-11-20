"""Correct the server-rejected Task12 schema without changing fitted inference."""

from collections import Counter
from pathlib import Path
import csv
import datetime
import gzip
import hashlib
import itertools
import json
import os
import shutil
import tempfile


ALLOWED_MEIOSES = {
    'siblings': 2, 'avuncular': 3, 'avuncular_once_removed': 4,
    'full_1st_cousin': 4, 'full_1st_cousin_once_removed': 5,
    'full_1st_cousin_twice_removed': 6, 'full_2nd_cousin': 6,
    'full_2nd_cousin_once_removed': 7, 'full_2nd_cousin_twice_removed': 8,
    'full_3rd_cousin': 8, 'full_3rd_cousin_once_removed': 9,
    'full_4th_cousin': 10,
}
ORIGINAL_MEIOSES = ALLOWED_MEIOSES | {
    'full_1st_cousin_thrice_removed': 7, 'avuncular_twice_removed': 5,
}
REPLACEMENTS = {
    'full_1st_cousin_thrice_removed': 'full_2nd_cousin_once_removed',
    'avuncular_twice_removed': 'full_1st_cousin_once_removed',
}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, data):
    with path.open('x') as stream:
        json.dump(data, stream, indent=2)
        stream.write('\n')


def main():
    parent = Path('/private/tmp/ghist-final-submissions-20261120/task12')
    mirror = Path(__file__).resolve().parent
    original = parent / 'ready_for_root/task12.final.csv'
    manifest = parent / 'ready_for_root/candidate_manifest.json'
    vcf = parent / 'relatedness4.final.vcf.gz'
    output = parent / 'task12_final_schema_corrected.csv'
    provenance_path = parent / 'task12_final_schema_correction_provenance.json'
    invalid_path = parent / 'submission_receipt_invalid_9785363.json'
    if any(path.exists() for path in [output, provenance_path, invalid_path]):
        raise FileExistsError('Schema correction outputs already exist')
    bound = {str(path): sha256(path) for path in [original, manifest, vcf, Path(__file__).resolve()]}
    expected_original = 'b94c2bef386ffcf062ad2fa1fe0c0593b2d2dadfc49efdaa916596897f718dac'
    candidate = json.loads(manifest.read_text())
    if bound[str(original)] != expected_original or candidate['candidate_sha256'] != expected_original:
        raise ValueError('Unexpected original candidate bytes')
    if bound[str(vcf)] != candidate['input_sha256'] or candidate['all_pairs'] != 43956 or candidate['sample_count'] != 297:
        raise ValueError('Original input/coverage binding mismatch')
    roster = None
    with gzip.open(vcf, 'rt') as stream:
        for line in stream:
            if line.startswith('#CHROM'):
                roster = line.rstrip('\n').split('\t')[9:]
                break
            if not line.startswith('#'):
                raise ValueError('VCF data encountered before the sample header')
    if roster is None or len(roster) != 297 or len(set(roster)) != 297:
        raise ValueError('Expected the original full 297-person header')
    with original.open(newline='') as stream:
        reader = csv.DictReader(stream, strict=True)
        if reader.fieldnames != ['indID1', 'indID2', 'relation']:
            raise ValueError('Original must have exactly three columns')
        rows = list(reader)
    pairs = [tuple(sorted([row['indID1'], row['indID2']])) for row in rows]
    expected_pairs = set(itertools.combinations(sorted(roster), 2))
    if len(rows) != 43956 or len(set(pairs)) != 43956 or set(pairs) != expected_pairs:
        raise ValueError('Original full-roster pair coverage mismatch')
    for old, new in REPLACEMENTS.items():
        equivalents = [name for name, degree in ALLOWED_MEIOSES.items() if degree == ORIGINAL_MEIOSES[old]]
        if equivalents != [new]:
            raise ValueError('Replacement is not the unique accepted same-meiosis label')
    corrected, changed = [], Counter()
    for row in rows:
        if set(row) != {'indID1', 'indID2', 'relation'}:
            raise ValueError('Malformed original row')
        degree = ORIGINAL_MEIOSES[row['relation']]
        label = REPLACEMENTS.get(row['relation'], row['relation'])
        if row['relation'] != label:
            changed[row['relation']] += 1
        if label not in ALLOWED_MEIOSES or ALLOWED_MEIOSES[label] != degree:
            raise ValueError('Corrected row violates server label or count rules')
        corrected.append(row | {'relation': label, 'meioses_count': str(degree)})
    if dict(changed) != {'full_1st_cousin_thrice_removed': 4275, 'avuncular_twice_removed': 3025}:
        raise ValueError('Unsupported-label counts differ from the rejection evidence')
    descriptor, name = tempfile.mkstemp(prefix='.task12-schema-', suffix='.csv', dir=parent)
    temporary = Path(name)
    try:
        with os.fdopen(descriptor, 'w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=['indID1', 'indID2', 'relation', 'meioses_count'])
            writer.writeheader()
            writer.writerows(corrected)
        with temporary.open(newline='') as stream:
            checked = list(csv.DictReader(stream, strict=True))
        if checked != corrected:
            raise ValueError('Corrected CSV round-trip mismatch')
        for before, after in zip(rows, checked):
            if any(before[key] != after[key] for key in ['indID1', 'indID2']):
                raise ValueError('Original pair row/order changed')
            if int(after['meioses_count']) != ORIGINAL_MEIOSES[before['relation']]:
                raise ValueError('Original implied meiosis count changed')
        if any(sha256(path) != digest for path, digest in bound.items()):
            raise ValueError('Input changed during schema correction')
        os.link(temporary, output)
    finally:
        temporary.unlink(missing_ok=True)
    receipt = {
        'status': 'INVALID', 'submission_id': '9785363', 'queue_id': '9619678',
        'entity_id': 'syn78780267', 'entity_version': 1,
        'entity_name': 'Task12 frozen relatedness baseline 2026-11-20',
        'email_id': '1a1fae902fbf3f6d', 'email_time_utc': '2026-11-20T21:42:16Z',
        'evidence_source': 'Root-reported actual server rejection email; retained separately from corrected candidate',
        'not_scored': True, 'original_candidate_sha256': expected_original,
        'server_allowed_labels': list(ALLOWED_MEIOSES),
        'unexpected_labels': {
            'full_1st_cousin_thrice_removed': {'rows': 4275, 'first_line': 5},
            'avuncular_twice_removed': {'rows': 3025, 'first_line': 22},
        },
    }
    write_json(invalid_path, receipt)
    provenance = {
        'state': 'SERVER_SCHEMA_CORRECTED_PENDING_ROOT_RETRY',
        'created_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'original_columns': ['indID1', 'indID2', 'relation'],
        'output_columns': ['indID1', 'indID2', 'relation', 'meioses_count'],
        'meioses_column': 'Appended deterministically from each original label before schema merging',
        'original_label_to_meioses': ORIGINAL_MEIOSES,
        'server_allowed_label_to_meioses': ALLOWED_MEIOSES,
        'replacement_mapping': REPLACEMENTS, 'replacement_counts': dict(changed),
        'mapping_interpretation': 'Unique accepted same-meiosis label; different pedigree relationships are merged, not biologically equivalent',
        'pairs': 43956, 'people': 297,
        'all_original_ids_pair_order_and_implied_meioses_preserved': True,
        'all_other_labels_preserved': True,
        'vcf_access': 'Original SHA256 verification and sample header only; no genotype parsing or model execution',
        'no_model_probability_refit_or_score_feedback': True,
        'input_sha256_before_after': bound,
        'output_path': str(output), 'output_sha256': sha256(output),
        'invalid_receipt_sha256': sha256(invalid_path),
        'scientific_qualification': False, 'submission_performed': False,
    }
    write_json(provenance_path, provenance)
    for path in [output, provenance_path, invalid_path]:
        target = mirror / path.name
        if target.exists():
            raise FileExistsError(target)
        shutil.copy2(path, target)
        if sha256(target) != sha256(path):
            raise ValueError('Mirror hash mismatch')
    print(json.dumps({'output': str(output), 'sha256': sha256(output),
                      'provenance_sha256': sha256(provenance_path),
                      'invalid_receipt_sha256': sha256(invalid_path),
                      'replacement_counts': dict(changed)}))


if __name__ == '__main__':
    main()
