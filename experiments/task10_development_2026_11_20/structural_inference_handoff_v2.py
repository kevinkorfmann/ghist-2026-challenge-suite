"""Inference-only VCF adapter; fixed map domain, original IDs, complete CSV pairs."""
from __future__ import annotations

import argparse
import csv
import gzip
import importlib.metadata
import itertools
import json
import os
import platform
import subprocess
import sys
import tempfile
import time
from collections import Counter
from pathlib import Path

import joblib
from ghist_submission.validate import validate_relatedness
from structural_confidence import attach_global_confidence
from structural_pair_candidate import infer, pair_features, pair_predictions

from scripts.develop_task10_detected_ibd import (
    genotype_roster,
    genome_length,
    read_haplotype_segments,
    read_map,
    save_json,
    sha256,
    verify_inputs,
)
from scripts.prepare_hapibd_vcf import prepare_vcf


def crop_to_map(source, output, bounds, minimum_markers):
    """Include both frozen map endpoints; retain every in-domain allele verbatim."""
    read, kept, outside = Counter(), Counter(), Counter()
    previous, finished, active, samples = {}, set(), None, None
    descriptor, name = tempfile.mkstemp(prefix=f'.{output.name}.', suffix='.gz', dir=output.parent)
    os.close(descriptor)
    temporary = Path(name)
    opener = gzip.open if source.suffix == '.gz' else open
    try:
        with opener(source, 'rt') as reader, gzip.open(temporary, 'wt', compresslevel=1) as writer:
            for line in reader:
                if line.startswith('#CHROM\t'):
                    if samples is not None:
                        raise ValueError('Duplicate sample header')
                    samples = line.rstrip('\r\n').split('\t')[9:]
                    if len(set(samples)) != len(samples) or any(not value or value != value.strip() for value in samples):
                        raise ValueError('Invalid or duplicate sample identifiers')
                if line.startswith('#'):
                    writer.write(line)
                    continue
                if samples is None:
                    raise ValueError('Missing VCF sample header')
                fields = line.rstrip('\r\n').split('\t')
                if len(fields) != 9 + len(samples):
                    raise ValueError('Inconsistent genotype columns')
                chromosome, position = fields[0], int(fields[1])
                if position < 1 or position <= previous.get(chromosome, 0):
                    raise ValueError('Expected strictly increasing positive VCF positions')
                previous[chromosome] = position
                if chromosome != active:
                    if chromosome in finished:
                        raise ValueError('Repeated disjoint chromosome block')
                    if active is not None:
                        finished.add(active)
                    active = chromosome
                read[chromosome] += 1
                if chromosome not in bounds or not bounds[chromosome][0] <= position <= bounds[chromosome][1]:
                    outside[chromosome] += 1
                    continue
                writer.write(line)
                kept[chromosome] += 1
        if samples is None or any(kept[chromosome] < minimum_markers for chromosome in bounds):
            raise ValueError('No or inadequate in-domain markers on a required chromosome')
        os.replace(temporary, output)
    finally:
        temporary.unlink(missing_ok=True)
    return {'samples': samples, 'read_per_chromosome': dict(read), 'retained_per_chromosome': dict(kept),
        'excluded_per_chromosome': dict(outside), 'records_read': sum(read.values()),
        'records_retained': sum(kept.values()), 'records_excluded': sum(outside.values()),
        'inclusive_map_bounds': bounds, 'alleles_modified': 0,
        'scope': 'Fixed training-map domain only; no genetic-map extrapolation, imputation, phase inference or MAF tuning'}


def write_candidate_csv(rows, samples, output):
    expected = set(itertools.combinations(sorted(samples), 2))
    observed = {tuple(sorted([row['indID1'], row['indID2']])) for row in rows}
    if len(set(samples)) != len(samples) or observed != expected or len(rows) != len(expected):
        raise ValueError('Candidate must contain every unique input-roster pair exactly once')
    descriptor, name = tempfile.mkstemp(prefix=f'.{output.name}.', suffix='.csv', dir=output.parent)
    os.close(descriptor)
    temporary = Path(name)
    try:
        with temporary.open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=['indID1', 'indID2', 'relation'])
            writer.writeheader()
            writer.writerows({key: row[key] for key in writer.fieldnames} for row in rows)
        report = validate_relatedness(temporary, 10, expected)
        os.replace(temporary, output)
    finally:
        temporary.unlink(missing_ok=True)
    return report


def eligible_markers(path):
    counts = Counter()
    with gzip.open(path, 'rt') as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            fields = line.split('\t')
            if ',' not in fields[4] and fields[4] != '.':
                counts[fields[0]] += 1
    return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-vcf', type=Path, required=True)
    parser.add_argument('--input-sha256', required=True)
    parser.add_argument('--expected-people', type=int, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--handoff-source', type=Path, required=True)
    parser.add_argument('--missing-policy', choices=['error', 'drop-site'], default='error')
    parser.add_argument('--missing-diagnostic-report', type=Path)
    args = parser.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):
        parser.error('Genotype inference requires Slurm')
    source = args.handoff_source.resolve()
    manifest = json.loads((source / 'source_manifest.json').read_text())
    verify_inputs(manifest)
    protocol_path = source / 'structural_inference_handoff_protocol.json'
    protocol = json.loads(protocol_path.read_text())
    frozen = Path(protocol['original_source'])
    qualification = Path(protocol['qualification_source'])
    origins = {name: str(Path(sys.modules[name].__file__).resolve()) for name in protocol['import_origins']}
    if origins != protocol['import_origins']:
        raise ValueError('Imported modules do not originate from frozen source snapshots')
    verify_inputs(json.loads((frozen / 'source_manifest.json').read_text()))
    verify_inputs(json.loads((qualification / 'source_manifest.json').read_text()))
    if sha256(args.input_vcf) != args.input_sha256:
        raise ValueError('Input VCF differs from root-frozen input SHA256')
    runtime = {'python': platform.python_version(), 'packages': {
        name: importlib.metadata.version(name) for name in protocol['runtime']['packages']}}
    if runtime != protocol['runtime']:
        raise ValueError('Inference runtime differs from frozen model runtime')
    if int(os.environ.get('SLURM_CPUS_PER_TASK', 0)) != protocol['cpus']:
        raise ValueError('Unexpected allocated CPU count')
    args.output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    confidence_path, model_path = map(Path, [protocol['confidence_path'], protocol['model_path']])
    freeze_path = Path(protocol['model_calibration_path'])
    bound = {str(path): sha256(path) for path in [args.input_vcf, confidence_path, model_path, freeze_path]}
    if args.missing_policy == 'drop-site':
        if not args.missing_diagnostic_report:
            raise ValueError('Optional complete-site policy requires passed frozen diagnostic')
        diagnostic = json.loads(args.missing_diagnostic_report.read_text())
        if diagnostic.get('all_gates_pass') is not True or diagnostic.get('protocol_sha256') != protocol['pragmatic_missingness_protocol_sha256']:
            raise ValueError('Optional missingness diagnostic has not passed the fixed policy')
        verify_inputs(diagnostic['input_sha256'])
        bound[str(args.missing_diagnostic_report)] = sha256(args.missing_diagnostic_report)
    confidence = json.loads(confidence_path.read_text())
    freeze = json.loads(freeze_path.read_text())
    if sha256(model_path) != protocol['model_sha256'] or sha256(model_path) != freeze['model_hashes']['pair_only']:
        raise ValueError('Pair model changed')
    if sha256(confidence_path) != protocol['confidence_sha256'] or confidence['confidence'] != protocol['global_correctness_confidence']:
        raise ValueError('Global confidence changed')
    save_json(args.output / 'initial_provenance.json', {'input_sha256': args.input_sha256,
        'source_manifest_sha256': sha256(source / 'source_manifest.json'),
        'protocol_sha256': sha256(protocol_path), 'runtime': runtime,
        'import_origins': origins, 'missing_policy': args.missing_policy,
        'slurm_job_id': os.environ['SLURM_JOB_ID'], 'source_readiness': 'Inference handoff only; root makes submission decision'})
    cropped, prepared = [args.output / name for name in ['map_domain.vcf.gz', 'prepared.vcf.gz']]
    crop = crop_to_map(args.input_vcf, cropped, protocol['map_bounds'], protocol['minimum_markers_per_chromosome'])
    save_json(args.output / 'map_domain_audit.json', crop)
    samples = crop['samples']
    if args.expected_people < 2 or len(samples) != args.expected_people:
        raise ValueError('Input roster differs from the complete sample count declared by root')
    eligible = eligible_markers(cropped)
    contract = prepare_vcf(cropped, prepared, missing_policy=args.missing_policy)
    save_json(args.output / 'phase_missingness_audit.json', contract)
    if contract['imputed_alleles'] or (args.missing_policy == 'error' and contract['missing_records_dropped']) or genotype_roster(prepared) != samples:
        raise ValueError('Input phase/missingness/roster contract changed')
    # Count retained biallelic records after the frozen preparer's own filtering.
    retained = Counter()
    with gzip.open(prepared, 'rt') as handle:
        for line in handle:
            if not line.startswith('#'):
                retained[line.split('\t', 1)[0]] += 1
    save_json(args.output / 'prepared_markers.json', dict(retained))
    if any(retained[chromosome] < protocol['minimum_markers_per_chromosome'] for chromosome in protocol['map_bounds']):
        raise ValueError('Inadequate retained biallelic markers after preparation')
    fractions = {chromosome: retained[chromosome] / eligible[chromosome] for chromosome in protocol['map_bounds']}
    save_json(args.output / 'retention_audit.json', {'eligible': dict(eligible), 'retained': dict(retained), 'fractions': fractions, 'missing_policy': args.missing_policy})
    if args.missing_policy == 'drop-site' and min(fractions.values()) < .75:
        raise ValueError('Optional complete-site filter retains less than75% on a chromosome')
    assets = frozen / 'assets'
    map_path = assets / 'external/ped-sim-map/refined_average.plink.map'
    if sha256(map_path) != protocol['map_sha256']:
        raise ValueError('Frozen genetic map changed')
    genetic_map = read_map(map_path)
    if {chrom: [int(pos[0]), int(pos[-1])] for chrom, (pos, _) in genetic_map.items()} != protocol['map_bounds']:
        raise ValueError('Crop bounds differ from exact model map')
    command = ['java', '-Xmx28g', '-jar', str(assets / 'external/hap-ibd/hap-ibd.jar'),
        f'gt={prepared}', f'map={map_path}', f'out={args.output / "detected"}',
        'min-seed=1.0', 'min-extend=0.2', 'min-output=2.0', 'min-markers=100', 'max-gap=1000',
        f'nthreads={protocol["cpus"]}']
    with (args.output / 'detector.log').open('w') as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=False)
    save_json(args.output / 'detector_command.json', {'command': command, 'exit_code': result.returncode})
    if result.returncode:
        raise RuntimeError('Frozen hap-IBD detector failed')
    grouped, _ = read_haplotype_segments(args.output / 'detected.ibd.gz')
    segments = {}
    for (a, b, chromosome), rows in grouped.items():
        if a not in samples or b not in samples or a == b:
            raise ValueError('Detector contains an unexpected or self pair')
        segments.setdefault((a, b), {})[chromosome] = rows
    features = pair_features(samples, segments, genetic_map, genome_length(map_path))
    save_json(args.output / 'features.json', features)
    bundle = joblib.load(model_path)
    pair_rows = pair_predictions(features, bundle, freeze['temperatures']['pair_only']['temperature'])
    save_json(args.output / 'pair_only_predictions.json', pair_rows)
    pair_report = write_candidate_csv(pair_rows, samples, args.output / 'pair_only.csv')
    structural, tree = infer(pair_rows)
    predictions = attach_global_confidence(structural, list(map(str, bundle['model'].classes_)), confidence)
    save_json(args.output / 'structural_predictions.json', predictions)
    save_json(args.output / 'structural_tree.json', tree)
    candidate = args.output / 'candidate.csv'
    report = write_candidate_csv(predictions, samples, candidate)
    verify_inputs(bound)
    verify_inputs(manifest)
    save_json(args.output / 'result.json', {'state': 'COMPLETE_INFERENCE_CANDIDATE',
        'people': len(samples), 'expected_people_from_root': args.expected_people, 'expected_pairs': args.expected_people * (args.expected_people - 1) // 2, 'pairs': len(predictions), 'csv_validation': report,
        'pair_only_csv_validation': pair_report, 'input_and_model_confidence_sha256': bound,
        'candidate_sha256': sha256(candidate), 'seconds_before_artifact_hashing': time.monotonic() - started,
        'protocol_sha256': sha256(protocol_path), 'submission_performed': False,
        'strict_qualification_status': 'Refer to independently scored synthetic qualification; this run does not change it',
        'confidence_interpretation': 'Constant global correctness estimate; zero individual resolution; not root posterior',
        'output_sha256': {str(path): sha256(path) for path in args.output.iterdir() if path.is_file()}})


if __name__ == '__main__':
    main()
