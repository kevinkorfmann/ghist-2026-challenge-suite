"""Fixed 25% site thinning on retained development genotypes; no fitting."""
from __future__ import annotations
import argparse
import gzip
import json
import os
import subprocess
from collections import Counter
from pathlib import Path
import joblib
import numpy as np
from structural_pair_candidate import infer, pair_features, pair_predictions
from scripts.develop_task10_detected_ibd import (genotype_roster, genome_length, read_haplotype_segments, read_map, save_json, sha256, verify_inputs)


def thin(source, destination, index):
    counts = Counter()
    with gzip.open(source, 'rt') as handle:
        for line in handle:
            if not line.startswith('#'):
                counts[line.split('\t', 1)[0]] += 1
    removed, seeds = {}, {}
    for chromosome, count in counts.items():
        seed = int(np.random.SeedSequence([2026207010, index, int(chromosome)]).generate_state(1)[0])
        seeds[chromosome] = seed
        removed[chromosome] = set(map(int, np.random.default_rng(seed).choice(count, count // 4, replace=False)))
    positions, kept = Counter(), Counter()
    with gzip.open(source, 'rt') as reader, gzip.open(destination, 'wt', compresslevel=1) as writer:
        for line in reader:
            if line.startswith('#'):
                writer.write(line)
                continue
            chromosome = line.split('\t', 1)[0]
            if positions[chromosome] not in removed[chromosome]:
                writer.write(line)
                kept[chromosome] += 1
            positions[chromosome] += 1
    return {'read': dict(counts), 'retained': dict(kept), 'seeds': seeds, 'alleles_modified': 0}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--index', type=int, required=True)
    parser.add_argument('--root', type=Path, required=True)
    args = parser.parse_args()
    if not os.environ.get('SLURM_JOB_ID') or int(os.environ.get('SLURM_CPUS_PER_TASK', '0')) != 4:
        parser.error('Requires scheduled four-CPU allocation')
    protocol = json.loads((args.source / 'pragmatic_missingness_protocol.json').read_text())
    manifest = json.loads((args.source / 'source_manifest.json').read_text())
    verify_inputs(manifest)
    if args.index not in protocol['development_cohorts']:
        raise ValueError('Undeclared development cohort')
    frozen = args.root / 'detected-panel-source'
    verify_inputs(json.loads((frozen / 'source_manifest.json').read_text()))
    original = args.root / 'detected_panel_v2' / f'cohort_{args.index:03d}'
    evaluation = args.root / 'detected_panel_v2/model_evaluation'
    inputs = {str(path): sha256(path) for path in [original / 'prepared.vcf.gz', evaluation / 'pair_only_model.joblib', evaluation / 'model_calibration_freeze.json']}
    if inputs[str(evaluation / 'pair_only_model.joblib')] != protocol['pair_model_sha256']:
        raise ValueError('Model changed')
    output = args.root / 'site_thinning_diagnostic_v1' / f'cohort_{args.index:03d}'
    output.mkdir(parents=True, exist_ok=False)
    prepared = output / 'thinned.vcf.gz'
    save_json(output / 'thinning.json', thin(original / 'prepared.vcf.gz', prepared, args.index))
    assets = frozen / 'assets'
    map_path = assets / 'external/ped-sim-map/refined_average.plink.map'
    command = ['java', '-Xmx14g', '-jar', str(assets / 'external/hap-ibd/hap-ibd.jar'), f'gt={prepared}', f'map={map_path}', f'out={output / "detected"}', 'min-seed=1.0', 'min-extend=0.2', 'min-output=2.0', 'min-markers=100', 'max-gap=1000', 'nthreads=4']
    with (output / 'detector.log').open('w') as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
    grouped, _ = read_haplotype_segments(output / 'detected.ibd.gz')
    segments = {}
    for (a,b,chromosome), rows in grouped.items():
        segments.setdefault((a,b), {})[chromosome] = rows
    samples = genotype_roster(prepared)
    features = pair_features(samples, segments, read_map(map_path), genome_length(map_path))
    bundle = joblib.load(evaluation / 'pair_only_model.joblib')
    freeze = json.loads((evaluation / 'model_calibration_freeze.json').read_text())
    pair_rows = pair_predictions(features, bundle, freeze['temperatures']['pair_only']['temperature'])
    rows, tree = infer(pair_rows)
    save_json(output / 'predictions.json', rows)
    save_json(output / 'tree.json', tree)
    save_json(output / 'pair_predictions.json', pair_rows)
    save_json(output / 'samples.json', samples)
    verify_inputs(inputs)
    verify_inputs(manifest)
    save_json(output / 'result.json', {'state': 'COMPLETE_DEVELOPMENT_DIAGNOSTIC', 'input_sha256': inputs, 'source_manifest_sha256': sha256(args.source / 'source_manifest.json'), 'protocol_sha256': sha256(args.source / 'pragmatic_missingness_protocol.json'), 'index': args.index, 'pairs': len(rows), 'slurm_job_id': os.environ['SLURM_JOB_ID'], 'output_sha256': {str(p): sha256(p) for p in output.iterdir() if p.is_file()}})


if __name__ == '__main__':
    main()
