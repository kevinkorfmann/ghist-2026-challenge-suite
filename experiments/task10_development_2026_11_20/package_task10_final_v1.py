"""Package frozen predictions only; no VCF access, inference, or model selection."""

from pathlib import Path
import csv
import datetime
import hashlib
import itertools
import json
import os
import sys
import tempfile


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    project = Path('/vast/projects/smathi/cohort/korfmann/ghist2026')
    final = project / 'final_runs/task10-root-20261120-2028'
    run = final / 'pragmatic_fallback_v1_run'
    frozen = project / 'experiments/task10-development-20261120/detected-panel-source'
    source_manifest = frozen / 'source_manifest.json'
    if sha256(source_manifest) != '7d4c00c9cfd9d197fce722e5de566078336ff0e6c2d7a39a84504e28a3bd9489':
        raise ValueError('Unexpected frozen source manifest')
    sources = json.loads(source_manifest.read_text())
    source_paths = [frozen / 'src' / name for name in [
        'ghist_relatedness/__init__.py', 'ghist_relatedness/labels.py',
        'ghist_submission/__init__.py', 'ghist_submission/validate.py']]
    for path in source_paths:
        if sha256(path) != sources[str(path)]:
            raise ValueError(f'Frozen import changed: {path}')
    sys.path.insert(0, str(frozen / 'src'))
    from ghist_relatedness import labels
    from ghist_submission import validate
    if Path(labels.__file__).resolve() != source_paths[1] or Path(validate.__file__).resolve() != source_paths[3]:
        raise ValueError('Unexpected frozen import origin')
    paths = [run / 'candidate.csv', run / 'structural_predictions.json',
             run / 'initial_provenance.json', run / 'result.json',
             final / 'vcf_preflight.json', source_manifest, Path(__file__).resolve(), *source_paths]
    bound = {str(path): sha256(path) for path in paths}
    result = json.loads((run / 'result.json').read_text())
    initial = json.loads((run / 'initial_provenance.json').read_text())
    preflight = json.loads((final / 'vcf_preflight.json').read_text())
    roster = preflight['sample_ids']
    if len(roster) != 263 or len(set(roster)) != 263 or preflight['sample_count'] != 263:
        raise ValueError('Full initial 263-person roster required')
    if preflight['input_sha256'] != initial['input_sha256']:
        raise ValueError('Preflight and inference inputs differ')
    if result['state'] != 'COMPLETE_INFERENCE_CANDIDATE' or result['people'] != 263 or result['pairs'] != 34453:
        raise ValueError('Inference did not complete the full roster')
    for path in [run / 'candidate.csv', run / 'structural_predictions.json', run / 'initial_provenance.json']:
        if bound[str(path)] != result['output_sha256'][str(path)]:
            raise ValueError('Final output binding mismatch')
    mapping = {label: up + down for (up, down, _), label in labels.HUMAN_LABELS.items()}
    if mapping != validate.relation_meioses(10) or len(mapping) != 14:
        raise ValueError('Frozen meiosis mappings disagree')
    expected = set(itertools.combinations(sorted(roster), 2))
    if len(expected) != 34453:
        raise ValueError('Unexpected pair denominator')
    original_validation = validate.validate_relatedness(run / 'candidate.csv', 10, expected)
    with (run / 'candidate.csv').open(newline='') as stream:
        reader = csv.DictReader(stream, strict=True)
        if reader.fieldnames != ['indID1', 'indID2', 'relation']:
            raise ValueError('Unexpected raw candidate header')
        rows = list(reader)
    predictions = json.loads((run / 'structural_predictions.json').read_text())
    projected = [{key: row[key] for key in ['indID1', 'indID2', 'relation']} for row in predictions]
    if rows != projected:
        raise ValueError('Raw CSV and frozen prediction rows differ')
    target = final / 'task10_final_candidate.csv'
    provenance_path = final / 'task10_final_packaging_provenance.json'
    if target.exists() or provenance_path.exists():
        raise FileExistsError('Packaging outputs already exist')
    descriptor, temporary_name = tempfile.mkstemp(prefix='.task10-package-', suffix='.csv', dir=final)
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, 'w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=['indID1', 'indID2', 'relation', 'meioses_count'])
            writer.writeheader()
            writer.writerows(row | {'meioses_count': mapping[row['relation']]} for row in rows)
        validation = validate.validate_relatedness(temporary, 10, expected)
        with temporary.open(newline='') as stream:
            packaged = list(csv.DictReader(stream, strict=True))
        if [{key: row[key] for key in ['indID1', 'indID2', 'relation']} for row in packaged] != rows:
            raise ValueError('Packaging altered row order, IDs, or relations')
        for path, digest in bound.items():
            if sha256(path) != digest:
                raise ValueError('Packaging input changed during validation')
        os.link(temporary, target)
    finally:
        temporary.unlink(missing_ok=True)
    provenance = {
        'state': 'PACKAGED_COMPLETE_FROZEN_PREDICTIONS',
        'created_at_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'scope': 'Append deterministic meioses_count only; no VCF access, inference, tuning or submission',
        'rows_preserved_exactly_in_order': True,
        'all_263_roster_pairs_preserved': True,
        'mapping_rule': 'Frozen HUMAN_LABELS tuple up + down; identical to relation_meioses(10)',
        'label_to_meioses_count': mapping,
        'original_validation': original_validation,
        'packaged_validation': validation,
        'input_sha256': bound,
        'output_path': str(target),
        'output_sha256': sha256(target),
        'strict_qualification_status': 'FAILED_UNCHANGED; packaging is not scientific qualification',
        'submission_performed': False,
    }
    with provenance_path.open('x') as stream:
        json.dump(provenance, stream, indent=2)
        stream.write('\n')
    print(json.dumps({'output': str(target), 'sha256': sha256(target), 'validation': validation,
                      'provenance': str(provenance_path), 'provenance_sha256': sha256(provenance_path)}))


if __name__ == '__main__':
    main()
