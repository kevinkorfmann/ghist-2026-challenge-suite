"""Separately declared practical scaling; unchanged candidate, no strict qualification claim."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import joblib
from detected_family_panel import select_observed
from nested_qualification_design import compile_definition
from structural_confidence import attach_global_confidence
from structural_pair_candidate import pair_features, pair_predictions
from structural_fallback_policy import infer_with_fallback as infer

from scripts.develop_task10_detected_ibd import (
    complete_features,
    derived_seed,
    detected_segments,
    family_truth,
    frozen_pair_predictions,
    generate_founders,
    genome_length,
    map_bounds,
    parents_from_fam,
    read_haplotype_segments,
    read_map,
    save_json,
    sha256,
    verify_inputs,
    write_features,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--frozen-source', type=Path, required=True)
    parser.add_argument('--qualification-source', type=Path, required=True)
    parser.add_argument('--launch-root', type=Path, required=True)
    parser.add_argument('--out-root', type=Path, required=True)
    parser.add_argument('--index', type=int, required=True)
    args = parser.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):
        parser.error('Qualification requires Slurm')
    if args.index not in [8, 9]:
        raise ValueError('Practical scaling permits only the same two fixed large cohorts')
    if datetime.now(timezone.utc) >= datetime(2026, 11, 20, 21, 50, tzinfo=timezone.utc):
        raise ValueError('Latest practical scaling start passed')
    frozen, qualification = args.frozen_source, args.qualification_source
    manifest = json.loads((qualification / 'source_manifest.json').read_text())
    verify_inputs(manifest)
    verify_inputs(json.loads((frozen / 'source_manifest.json').read_text()))
    practical_source = Path(__file__).resolve().parent
    verify_inputs(json.loads((practical_source / 'source_manifest.json').read_text()))
    launch = json.loads((args.launch_root / 'launch_gate.json').read_text())
    if sha256(practical_source / 'fallback_scaling_protocol.json') != launch['practical_protocol_sha256']:
        raise ValueError('Practical protocol changed')
    verify_inputs(launch['input_sha256'])
    if not all(launch.get(key) is True for key in ['structural_development_gates_pass',
            'independent_development_rescore_verified', 'direct_feature_and_model_equivalence_pass',
            'runtime_preflight_pass', 'global_confidence_rule_frozen', 'source_and_plan_reviewed']):
        raise ValueError('Qualification prerequisites have not all passed')
    confidence_path = args.launch_root / 'confidence_freeze.json'
    if sha256(confidence_path) != launch['confidence_freeze_sha256']:
        raise ValueError('Confidence changed after launch freeze')
    confidence = json.loads(confidence_path.read_text())
    protocol_path = qualification / 'structural_qualification_protocol_v3.json'
    protocol = json.loads(protocol_path.read_text())
    if sha256(protocol_path) != confidence['protocol_sha256'] or sha256(protocol_path) != launch['protocol_sha256']:
        raise ValueError('Qualification policy differs from calibration/launch freeze')
    if sha256(frozen / 'source_manifest.json') != protocol['candidate_source_manifest_sha256']:
        raise ValueError('Original candidate source changed')
    item = protocol['cohorts'][args.index]
    if item['index'] != args.index:
        raise ValueError('Wrong qualification cohort')
    if launch.get('practical_small_all_gates_pass') is not True or launch.get('independent_small_raw_audit_pass') is not True:
        raise ValueError('Practical scaling requires fixed practical gates and independent small audit')
    if launch.get('strict_small_result_preserved') is not True:
        raise ValueError('Original strict failure or success must remain explicit')
    workers = int(os.environ.get('SLURM_CPUS_PER_TASK', '1'))
    if workers != (8 if args.index >= 8 else 4):
        raise ValueError('Resource allocation differs from prospective plan')
    definition, roster, expected = compile_definition(item['terminal_paths'], item['name'], item['people'])
    if definition != item['definition'] or roster != item['declared_roster']:
        raise ValueError('Qualification definition differs from prospective design')
    seed = item['seed']
    if derived_seed(seed, 0, 3) != item['derived_seeds']['meioses']:
        raise ValueError('Meiosis stream changed')
    for chrom, streams in item['derived_seeds']['chromosomes'].items():
        for stage, name in [(1, 'ancestry'), (2, 'mutation'), (4, 'dropout_contract_only')]:
            if derived_seed(seed, int(chrom), stage) != streams[name]:
                raise ValueError('Chromosome stream changed')
    output = args.out_root / f'cohort_{args.index:03d}'
    output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    save_json(output / 'preregistered_design.json', {'item': item, 'protocol_sha256': sha256(protocol_path),
        'launch_gate': launch, 'source_manifest_sha256': sha256(qualification / 'source_manifest.json')})
    commands = []

    def run(name, command):
        words = list(map(str, command))
        with (output / f'{name}.log').open('w') as log:
            result = subprocess.run(words, stdout=log, stderr=subprocess.STDOUT, check=False)
        commands.append({'name': name, 'command': words, 'exit_code': result.returncode})
        save_json(output / 'commands.json', commands)
        if result.returncode:
            raise RuntimeError(f'{name} failed')

    assets = frozen / 'assets'
    map_path = assets / 'external/ped-sim-map/refined_average.plink.map'
    definition_path = output / 'family.def'
    definition_path.write_text(definition)
    common = [assets / 'external/ped-sim/ped-sim', '-d', definition_path,
        '-m', assets / 'external/ped-sim-map/refined_mf.simmap',
        '--intf', assets / 'external/ped-sim/interfere/nu_p_campbell.tsv',
        '--seed', derived_seed(seed, 0, 3)]
    run('structure', [*common, '-o', output / 'structure', '--dry_run'])
    pedigree = parents_from_fam(output / 'structure-everyone.fam')
    founder_count = sum(not parents for parents in pedigree.values())
    founders = generate_founders(map_bounds(map_path), seed, founder_count, output / 'founders', workers)
    run('pedsim', [*common, '-i', founders, '-o', output / 'family', '--fam', '--seg_haps',
        '--keep_phase', '--err_rate', 0, '--miss_rate', 0])
    parents = parents_from_fam(output / 'family-everyone.fam')
    normalized = output / 'family.normalized.vcf.gz'
    identifiers = select_observed(output / 'family.vcf.gz', normalized, expected, parents)
    truth = family_truth(parents, identifiers)
    declared = {(identifiers[a], identifiers[b]): label for (a, b), label in expected.items()}
    if {(row['indID1'], row['indID2']): row['relation'] for row in truth} != declared:
        raise ValueError('Simulated parentage does not match predeclared truth')
    save_json(output / 'truth.json', truth)
    save_json(output / 'identifiers.json', identifiers)
    prepared = output / 'prepared.vcf.gz'
    run('prepare', [sys.executable, frozen / 'scripts/prepare_hapibd_vcf.py',
        '--input', normalized, '--output', prepared])
    run('detector', ['java', '-Xmx28g' if args.index >= 8 else '-Xmx14g', '-jar', assets / 'external/hap-ibd/hap-ibd.jar',
        f'gt={prepared}', f'map={map_path}', f'out={output / "detected"}', 'min-seed=1.0',
        'min-extend=0.2', 'min-output=2.0', 'min-markers=100', 'max-gap=1000', f'nthreads={workers}'])
    grouped, _ = read_haplotype_segments(output / 'detected.ibd.gz')
    segments = {}
    for (a, b, chrom), rows in grouped.items():
        segments.setdefault((a, b), {})[chrom] = rows
    samples = sorted(identifiers.values())
    features = pair_features(samples, segments, read_map(map_path), genome_length(map_path))
    save_json(output / 'features.json', features)
    evaluation = frozen.parent / 'detected_panel_v2/model_evaluation'
    freeze = json.loads((evaluation / 'model_calibration_freeze.json').read_text())
    model_path = evaluation / 'pair_only_model.joblib'
    if sha256(model_path) != freeze['model_hashes']['pair_only']:
        raise ValueError('Model changed after calibration freeze')
    bundle = joblib.load(model_path)
    predictions = pair_predictions(features, bundle, freeze['temperatures']['pair_only']['temperature'])
    save_json(output / 'pair_only_predictions.json', predictions)
    structural, tree = infer(predictions)
    structural = attach_global_confidence(structural, list(map(str, bundle['model'].classes_)), confidence)
    save_json(output / 'structural_predictions.json', structural)
    save_json(output / 'structural_tree.json', tree)
    raw = complete_features(detected_segments(output / 'detected.ibd.gz', map_path, 'raw_haplotype'),
        samples, genome_length(map_path))
    baseline_features = output / 'baseline_features.tsv'
    write_features(raw, baseline_features)
    frozen_pair_predictions(baseline_features,
        assets / 'artifacts/task10_large_family/clean_spectrum_v1/model.joblib',
        assets / 'artifacts/task10_large_family/pair_model_min5/model.joblib', output / 'baseline_predictions.json')
    run('missingness_contract', [sys.executable, qualification / 'verify_missingness_contract.py',
        '--protocol', protocol_path, '--cohort', output, '--index', args.index,
        '--output', output / 'missingness_contract'])
    verify_inputs(manifest)
    verify_inputs(launch['input_sha256'])
    save_json(output / 'result.json', {'state': 'COMPLETE_NEW_FALLBACK_POLICY_PROSPECTIVE_SCALING_DATA',
        'item': item, 'people': len(samples), 'pairs': len(features), 'founder_count': founder_count,
        'total_seconds_before_artifact_hashing': time.monotonic() - started, 'submission_ready': False, 'strict_qualification_claim': False,
        'slurm_job_id': os.environ['SLURM_JOB_ID'], 'cpus': workers,
        'slurm_array_job_id': os.environ.get('SLURM_ARRAY_JOB_ID'),
        'slurm_array_task_id': os.environ.get('SLURM_ARRAY_TASK_ID'),
        'runtime_memory_gate_requires_final_slurm_accounting': True,
        'missingness_contract_check': 'Verified explicit rejection only; no biological missingness robustness',
        'output_sha256': {str(path): sha256(path) for path in output.rglob('*') if path.is_file()}})


if __name__ == '__main__':
    main()
