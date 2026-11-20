"""Unlaunched qualification runner; requires separately frozen sources and models."""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import joblib
from calibrate_phase import temperature
from detected_family_panel import select_observed
from fast_multiway_features import family_features
from nested_qualification_design import compile_definition

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
    parser.add_argument('--out-root', type=Path, required=True)
    parser.add_argument('--index', type=int, required=True)
    args = parser.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):
        parser.error('Qualification requires Slurm')
    frozen, qualification = args.frozen_source, args.qualification_source
    manifest = json.loads((qualification / 'source_manifest.json').read_text())
    verify_inputs(manifest)
    verify_inputs(json.loads((frozen / 'source_manifest.json').read_text()))
    launch = json.loads((qualification / 'launch_gate.json').read_text())
    if not all(launch.get(key) is True for key in ['v2_development_gates_pass',
            'independent_v2_rescore_verified', 'feature_equivalence_pass', 'model_equivalence_pass',
            'runtime_preflight_pass']):
        raise ValueError('Qualification prerequisites have not all passed')
    protocol_path = qualification / 'nested_qualification_protocol.json'
    protocol = json.loads(protocol_path.read_text())
    if sha256(frozen / 'source_manifest.json') != protocol['candidate_source_manifest_sha256']:
        raise ValueError('Original candidate source changed')
    item = protocol['cohorts'][args.index]
    if item['index'] != args.index:
        raise ValueError('Wrong qualification cohort')
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
    founders = generate_founders(map_bounds(map_path), seed, founder_count, output / 'founders', 4)
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
    run('detector', ['java', '-Xmx14g', '-jar', assets / 'external/hap-ibd/hap-ibd.jar',
        f'gt={prepared}', f'map={map_path}', f'out={output / "detected"}', 'min-seed=1.0',
        'min-extend=0.2', 'min-output=2.0', 'min-markers=100', 'max-gap=1000', 'nthreads=4'])
    grouped, _ = read_haplotype_segments(output / 'detected.ibd.gz')
    segments = {}
    for (a, b, chrom), rows in grouped.items():
        segments.setdefault((a, b), {})[chrom] = rows
    samples = sorted(identifiers.values())
    features, groups, _ = family_features(samples, segments, read_map(map_path), genome_length(map_path),
        joblib.load(frozen / 'models/phase_model.joblib'))
    save_json(output / 'features.json', features)
    save_json(output / 'inferred_groups.json', groups)
    freeze = json.loads((qualification / 'models/model_calibration_freeze.json').read_text())
    for variant in ['context', 'pair_only']:
        model_path = qualification / f'models/{variant}_model.joblib'
        if sha256(model_path) != freeze['model_hashes'][variant]:
            raise ValueError('Model changed after calibration freeze')
        bundle = joblib.load(model_path)
        values = bundle['model'].predict_proba([[row[name] for name in bundle['feature_names']] for row in features])
        values = temperature(values, freeze['temperatures'][variant]['temperature'])
        predictions = [{key: row[key] for key in ['indID1', 'indID2']} |
            {'relation': str(bundle['model'].classes_[probabilities.argmax()]),
                'probabilities': dict(zip(map(str, bundle['model'].classes_), map(float, probabilities), strict=True))}
            for row, probabilities in zip(features, values, strict=True)]
        save_json(output / f'{variant}_predictions.json', predictions)
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
    save_json(output / 'result.json', {'state': 'COMPLETE_CLEAN_GENOTYPE_QUALIFICATION_DATA',
        'item': item, 'people': len(samples), 'pairs': len(features), 'founder_count': founder_count,
        'total_seconds': time.monotonic() - started, 'submission_ready': False,
        'missingness_contract_check': 'Verified explicit rejection only; no biological missingness robustness',
        'output_sha256': {str(path): sha256(path) for path in output.rglob('*') if path.is_file()}})


if __name__ == '__main__':
    main()
