"""Fixed historical Task12 final inference; never trains or submits."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import platform
import re
import subprocess
import sys
import tempfile
from datetime import datetime, timezone

BASE = Path(__file__).resolve().parent
sys.path[:0] = [str(BASE), str(BASE / 'src')]


def digest(path, algorithm='sha256'):
    h = hashlib.new(algorithm)
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def write_json(path, payload):
    path = Path(path)
    fd, name = tempfile.mkstemp(prefix=path.name + '.', dir=path.parent)
    try:
        with os.fdopen(fd, 'w') as handle:
            json.dump(payload, handle, indent=2, allow_nan=False)
            handle.write('\n')
        os.replace(name, path)
    finally:
        Path(name).unlink(missing_ok=True)


def utc():
    return datetime.now(timezone.utc).isoformat()


def read_bound_json(path):
    data = Path(path).read_bytes()
    return json.loads(data), hashlib.sha256(data).hexdigest()


def verify_static(lock, lock_hash, manifest, manifest_hash):
    if digest(BASE / 'method_lock.json') != lock_hash or digest(BASE / 'source_manifest.json') != manifest_hash:
        raise ValueError('Bound protocol/source manifest changed during execution')
    for rel, expected in manifest['files'].items():
        if digest(BASE / rel) != expected:
            raise ValueError(f'Source changed: {rel}')
    for name, item in lock['assets'].items():
        if digest(item['path']) != item['sha256']:
            raise ValueError(f'Asset mismatch: {name}')
    if platform.python_version() != lock['python']:
        raise ValueError('Python version mismatch')
    java = subprocess.run(['java', '-version'], capture_output=True, text=True, check=True).stderr
    if java != lock['java_version']:
        raise ValueError('Java runtime mismatch')
    for package, expected in lock['packages'].items():
        if importlib.metadata.version(package) != expected:
            raise ValueError(f'Runtime version mismatch: {package}')


def verify_input(path, binding, binding_hash, expected_entity):
    if digest(binding) != binding_hash:
        raise ValueError('Input binding changed during execution')
    info = json.loads(Path(binding).read_text())
    if info['entity_id'] != expected_entity or not isinstance(info['version'], int) or info['version'] < 1:
        raise ValueError('Input identity/version mismatch')
    if Path(path).stat().st_size != info['bytes'] or digest(path) != info['sha256'] or digest(path, 'md5') != info['md5']:
        raise ValueError('Input bytes/hash mismatch')
    return info


def inspect_vcf(path):
    samples = None
    contigs = {}
    seen = set()
    last_chrom = None
    last_pos = 0
    count = 0
    observed = {}
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as handle:
        for line in handle:
            if line.startswith('##contig=<'):
                name = re.search(r'(?:<|,)ID=([^,>]+)', line)
                length = re.search(r'(?:<|,)length=(\d+)', line, re.IGNORECASE)
                if name:
                    if name[1] in contigs:
                        raise ValueError('Duplicate contig header')
                    contigs[name[1]] = int(length[1]) if length else None
            elif line.startswith('#CHROM\t'):
                if samples is not None:
                    raise ValueError('Duplicate sample header')
                samples = line.rstrip('\r\n').split('\t')[9:]
                if len(samples) < 2 or len(samples) != len(set(samples)) or any(not s or s != s.strip() for s in samples):
                    raise ValueError('Expected at least two unique nonempty sample IDs')
            elif not line.startswith('#'):
                fields = line.rstrip('\r\n').split('\t')
                if samples is None or len(fields) != len(samples) + 9:
                    raise ValueError('Invalid VCF record/header')
                chrom, pos = fields[0], int(fields[1])
                if pos < 1 or (contigs.get(chrom) is not None and pos > contigs[chrom]):
                    raise ValueError('VCF position outside declared coordinate domain')
                if chrom != last_chrom:
                    if chrom in seen:
                        raise ValueError('Chromosome blocks are not contiguous')
                    seen.add(chrom)
                    last_chrom, last_pos = chrom, 0
                if pos <= last_pos:
                    raise ValueError('VCF positions must be unique and strictly increasing per chromosome')
                last_pos = pos
                count += 1
                observed[chrom] = {'last_position': pos, 'records': observed.get(chrom, {}).get('records', 0) + 1}
    if not count:
        raise ValueError('Empty VCF')
    return samples, {'records': count, 'samples': len(samples), 'contigs': contigs, 'observed': observed}


def inspect_map(path, vcf_contigs, genome_cm):
    ranges = {}
    previous = {}
    with Path(path).open() as handle:
        for line in handle:
            chrom, _, cm, bp = line.split()
            cm, bp = float(cm), int(bp)
            if not math.isfinite(cm) or cm < 0 or bp < 0:
                raise ValueError('Invalid map coordinates')
            if chrom in previous and (bp <= previous[chrom][0] or cm < previous[chrom][1]):
                raise ValueError('Map must increase by physical position and not decrease by genetic position')
            previous[chrom] = (bp, cm)
            bounds = ranges.setdefault(chrom, [cm, cm])
            bounds[0], bounds[1] = min(bounds[0], cm), max(bounds[1], cm)
    if not set(vcf_contigs) <= set(ranges):
        raise ValueError('VCF chromosome absent from genetic map')
    total = sum(high - low for low, high in ranges.values())
    if not math.isclose(total, genome_cm, rel_tol=1e-12, abs_tol=1e-8):
        raise ValueError(f'Map genome length differs from locked feature contract: {total}')
    return {'genome_cm': total, 'chromosomes': sorted(ranges)}


def read_segments(path, samples, coordinate_info):
    grouped = defaultdict(list)
    n = 0
    with gzip.open(path, 'rt') as handle:
        for line in handle:
            fields = line.split()
            if len(fields) != 8:
                raise ValueError('hap-IBD output must have eight columns')
            left, h1, right, h2, chrom, start, end, length = fields
            start, end, length = int(start), int(end), float(length)
            if left not in samples or right not in samples or left == right or h1 not in {'1', '2'} or h2 not in {'1', '2'}:
                raise ValueError('Invalid haplotype/sample pair')
            if chrom not in coordinate_info['observed'] or not 1 <= start <= end <= coordinate_info['observed'][chrom]['last_position'] or not math.isfinite(length) or length <= 0:
                raise ValueError('Invalid detected segment coordinates/length')
            grouped[tuple(sorted((left, right)))].append((chrom, length))
            n += 1
    return grouped, n


def feature_records(grouped, samples, reference, lock):
    from ghist_submission.validate import all_pairs
    from scripts.fit_human_background_corrected import correct_records
    from scripts.fit_human_family_relationships import summarize_raw_lengths
    present = [{'indID1': a, 'indID2': b, **summarize_raw_lengths(grouped[(a, b)], lock['genome_cm'])} for a, b in sorted(grouped)]
    corrected = correct_records(present, reference, lock['minimum_cm'], lock['genome_cm']) if present else []
    missing = sorted(all_pairs(samples) - set(grouped))
    corrected += [{'indID1': a, 'indID2': b, **summarize_raw_lengths([], lock['genome_cm'])} for a, b in missing]
    return sorted(corrected, key=lambda row: (row['indID1'], row['indID2'])), len(missing)


def check_predictions(rows, samples):
    from ghist_submission.validate import all_pairs
    from scripts.fit_human_family_relationships import LABELS
    expected = all_pairs(samples)
    seen = set()
    for row in rows:
        pair = (row['indID1'], row['indID2'])
        if pair not in expected or pair in seen or row['relation'] not in LABELS:
            raise ValueError('Prediction pair/label domain error')
        seen.add(pair)
        p = row['probabilities']
        if set(p) != set(LABELS) or any(not math.isfinite(x) or x < 0 or x > 1 for x in p.values()):
            raise ValueError('Invalid probabilities')
        if not math.isclose(sum(p.values()), 1, rel_tol=1e-10, abs_tol=1e-10) or p[row['relation']] != max(p.values()) or row['confidence'] != p[row['relation']]:
            raise ValueError('Prediction probability/maximum mismatch')
    if seen != expected:
        raise ValueError('Incomplete all-pair prediction coverage')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf', type=Path, required=True)
    parser.add_argument('--binding', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    lock, lock_hash = read_bound_json(BASE / 'method_lock.json')
    manifest, manifest_hash = read_bound_json(BASE / 'source_manifest.json')
    _, binding_hash = read_bound_json(args.binding)
    verify_static(lock, lock_hash, manifest, manifest_hash)
    identity = verify_input(args.vcf, args.binding, binding_hash, lock['input_entity'])
    args.out.mkdir(parents=True, exist_ok=False)
    state = {'status': 'RUNNING', 'started_utc': utc(), 'lock_sha256': lock_hash, 'source_manifest_sha256': manifest_hash, 'input_binding_sha256': binding_hash, 'input': identity, 'scientific_limitations': lock['scientific_limitations']}
    write_json(args.out / 'status.json', state)
    try:
        import joblib
        import numpy as np
        from ghist_submission.validate import all_pairs, validate_submission
        from scripts.prepare_hapibd_vcf import prepare_vcf
        from scripts.fit_human_background_corrected import corrected_feature_names, observed_predictions
        from scripts.fit_human_family_relationships import LABELS
        from scripts.validate_human_background_correction import read_pedsim_with_lengths
        samples, vcf_info = inspect_vcf(args.vcf)
        map_info = inspect_map(lock['assets']['map']['path'], vcf_info['observed'], lock['genome_cm'])
        write_json(args.out / 'input_preflight.json', {'vcf': vcf_info, 'sample_ids': samples, 'map': map_info})
        prepared = args.out / 'prepared.vcf.gz'
        preparation = prepare_vcf(args.vcf, prepared, missing_policy=lock['preparation']['missing_policy'])
        if preparation['records_written'] == 0:
            raise ValueError('No retained biallelic polymorphic sites')
        write_json(args.out / 'preparation.json', preparation)
        prefix = args.out / 'segments'
        cmd = ['java', '-Xmx32g', '-jar', lock['assets']['jar']['path'], f'gt={prepared}', f'map={lock["assets"]["map"]["path"]}', f'out={prefix}']
        cmd += [f'{key}={value}' for key, value in lock['detector'].items()]
        write_json(args.out / 'detector_command.json', {'argv': cmd})
        with (args.out / 'detector.log').open('w') as handle:
            subprocess.run(cmd, stdout=handle, stderr=subprocess.STDOUT, check=True)
        grouped, segment_count = read_segments(Path(str(prefix) + '.ibd.gz'), set(samples), vcf_info)
        reference, _ = read_pedsim_with_lengths(Path(lock['assets']['reference']['path']), lock['genome_cm'])
        records, missing = feature_records(grouped, samples, reference, lock)
        model = joblib.load(lock['assets']['model']['path'])
        baseline = joblib.load(lock['assets']['baseline']['path'])
        if model['feature_names'] != corrected_feature_names(lock['minimum_cm']) or model['corrected_model_weight'] != lock['corrected_model_weight']:
            raise ValueError('Historical model feature/weight contract mismatch')
        if set(model['corrected_payload']['model'].classes_) != set(LABELS) or set(baseline['model'].classes_) != set(LABELS):
            raise ValueError('Model class domain mismatch')
        for row in records:
            if any(not np.isfinite(float(v)) or float(v) < 0 for k, v in row.items() if k not in {'indID1', 'indID2'}):
                raise ValueError('Nonfinite/negative feature')
        rows = observed_predictions(model['corrected_payload'], baseline['model'], baseline['feature_names'], records, model['feature_names'], model['corrected_model_weight'])
        check_predictions(rows, samples)
        candidate = args.out / 'candidate.pending.csv'
        with candidate.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=['indID1', 'indID2', 'relation'], extrasaction='ignore')
            writer.writeheader()
            writer.writerows(rows)
        validation = validate_submission(candidate, 12, 'final', all_pairs(samples))
        write_json(args.out / 'predictions_with_uncalibrated_probabilities.json', rows)
        verify_input(args.vcf, args.binding, binding_hash, lock['input_entity'])
        verify_static(lock, lock_hash, manifest, manifest_hash)
        final = args.out / 'task12.final.csv'
        candidate.rename(final)
        state.update(status='COMPLETE', completed_utc=utc(), output_sha256=digest(final), output_bytes=final.stat().st_size, validation=validation, segment_count=segment_count, segment_present_pairs=len(grouped), zero_segment_pairs=missing, all_pairs=len(records), sample_count=len(samples), relation_counts=dict(Counter(row['relation'] for row in rows)), zero_segment_policy=lock['zero_segment_policy'])
        write_json(args.out / 'status.json', state)
    except BaseException as error:
        state.update(status='FAILED', completed_utc=utc(), error=repr(error))
        write_json(args.out / 'status.json', state)
        raise


if __name__ == '__main__':
    main()
