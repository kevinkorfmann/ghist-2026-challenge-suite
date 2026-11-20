"""Fixed detected-IBD family panel with topology-disjoint training/calibration/development."""
from __future__ import annotations

import argparse
import gzip
import itertools
import json
import os
import subprocess
import sys
from pathlib import Path

import joblib
import numpy as np
from calibrate_phase import metrics, temperature
from multiway_candidate import (
    family_features,
    score_family_rows,
    whole_family_definition,
)

from scripts.develop_task10_detected_ibd import (
    complete_features,
    derived_seed,
    detected_segments,
    family_truth,
    file_hashes,
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


def bound_inputs(args):
    paths = [args.protocol.resolve(), Path(__file__).resolve(), args.phase_model.resolve(),
        args.asset_root / 'external/ped-sim/ped-sim',
        args.asset_root / 'external/ped-sim-map/refined_mf.simmap',
        args.asset_root / 'external/ped-sim-map/refined_average.plink.map',
        args.asset_root / 'external/ped-sim/interfere/nu_p_campbell.tsv',
        args.asset_root / 'external/hap-ibd/hap-ibd.jar',
        args.asset_root / 'artifacts/task10_large_family/clean_spectrum_v1/model.joblib',
        args.asset_root / 'artifacts/task10_large_family/pair_model_min5/model.joblib']
    hashes = file_hashes(args.source_root)
    paths.extend(sorted(Path(__file__).parent.glob('*.py')))
    hashes.update({str(path.resolve()): sha256(path) for path in paths})
    return hashes


def select_observed(source, destination, expected, parents):
    wanted = sorted({sample for pair in expected for sample in pair})
    identifiers = {sample: f'sample_{index:03d}' for index, sample in enumerate(wanted)}
    if any(len(parents[sample]) != 2 for sample in wanted):
        raise ValueError('Observed roster contains unexpected founders')
    variants, indices = 0, None
    with gzip.open(source, 'rt') as reader, gzip.open(destination, 'wt', compresslevel=1) as writer:
        for line in reader:
            if line.startswith('#CHROM'):
                fields = line.rstrip().split('\t')
                if not set(wanted) <= set(fields[9:]):
                    raise ValueError('Missing declared genotype sample')
                indices = [fields.index(sample) for sample in wanted]
                writer.write('\t'.join(fields[:9] + list(identifiers.values())) + '\n')
            elif line.startswith('#'):
                writer.write(line)
            else:
                if indices is None:
                    raise ValueError('Missing genotype header')
                fields = line.rstrip().split('\t')
                calls = [fields[index] for index in indices]
                if any(call.split(':')[0] not in {'0|0', '0|1', '1|0', '1|1'} for call in calls):
                    raise ValueError('Unexpected genotype missingness or phase')
                writer.write('\t'.join(fields[:9] + calls) + '\n')
                variants += 1
    if not variants:
        raise ValueError('Empty genotype panel')
    return identifiers


def generate(args, protocol):
    item = protocol['cohorts'][args.index]
    if item['index'] != args.index:
        raise ValueError('Index does not match frozen roster')
    output = args.out_root / f"cohort_{args.index:03d}"
    output.mkdir(parents=True, exist_ok=False)
    hashes = bound_inputs(args)
    save_json(output / 'protocol.json', {'item': item, 'hashes': hashes, 'stage': 'DETECTED_PANEL_DEVELOPMENT'})
    commands = []

    def run(name, words):
        command = list(map(str, words))
        with (output / f'{name}.log').open('w') as log:
            result = subprocess.run(command, check=False, stdout=log, stderr=subprocess.STDOUT)
        commands.append({'name': name, 'command': command, 'exit_code': result.returncode})
        save_json(output / 'commands.json', commands)
        if result.returncode:
            raise RuntimeError(f'{name} failed')

    assets = args.asset_root
    map_path = assets / 'external/ped-sim-map/refined_average.plink.map'
    depths, seed = item['depths'], item['seed']
    founders = generate_founders(map_bounds(map_path), seed,
        2 + len(depths) * (max(depths) - 1), output / 'founders', 4)
    definition, _, expected = whole_family_definition([depths], 'panel')
    definition_path = output / 'family.def'
    definition_path.write_text(definition)
    run('pedsim', [assets / 'external/ped-sim/ped-sim', '-d', definition_path,
        '-m', assets / 'external/ped-sim-map/refined_mf.simmap', '-i', founders,
        '-o', output / 'family', '--intf', assets / 'external/ped-sim/interfere/nu_p_campbell.tsv',
        '--seed', derived_seed(seed, 0, 3), '--fam', '--seg_haps', '--keep_phase', '--err_rate', 0, '--miss_rate', 0])
    parents = parents_from_fam(output / 'family-everyone.fam')
    normalized = output / 'family.normalized.vcf.gz'
    identifiers = select_observed(output / 'family.vcf.gz', normalized, expected, parents)
    truth = family_truth(parents, identifiers)
    declared = {(identifiers[a], identifiers[b]): label for (a, b), label in expected.items()}
    if {(row['indID1'], row['indID2']): row['relation'] for row in truth} != declared:
        raise ValueError('Actual simulated parentage differs from frozen topology')
    save_json(output / 'identifiers.json', identifiers)
    save_json(output / 'truth.json', truth)
    prepared = output / 'prepared.vcf.gz'
    run('prepare', [sys.executable, args.source_root / 'scripts/prepare_hapibd_vcf.py',
        '--input', normalized, '--output', prepared])
    run('detector', ['java', '-Xmx14g', '-jar', assets / 'external/hap-ibd/hap-ibd.jar',
        f'gt={prepared}', f'map={map_path}', f'out={output / "detected"}', 'min-seed=1.0',
        'min-extend=0.2', 'min-output=2.0', 'min-markers=100', 'max-gap=1000', 'nthreads=4'])
    grouped, _ = read_haplotype_segments(output / 'detected.ibd.gz')
    segments = {}
    for (a, b, chrom), rows in grouped.items():
        segments.setdefault((a, b), {})[chrom] = rows
    roster = sorted(identifiers.values())
    features, groups, _ = family_features(roster, segments, read_map(map_path), genome_length(map_path), joblib.load(args.phase_model))
    save_json(output / 'features.json', features)
    save_json(output / 'inferred_groups.json', groups)
    raw = complete_features(detected_segments(output / 'detected.ibd.gz', map_path, 'raw_haplotype'), roster, genome_length(map_path))
    feature_path = output / 'baseline_features.tsv'
    write_features(raw, feature_path)
    frozen_pair_predictions(feature_path,
        assets / 'artifacts/task10_large_family/clean_spectrum_v1/model.joblib',
        assets / 'artifacts/task10_large_family/pair_model_min5/model.joblib',
        output / 'baseline_predictions.json')
    verify_inputs(hashes)
    retained = {str(path.resolve()): sha256(path) for path in output.rglob('*') if path.is_file()}
    save_json(output / 'result.json', {'state': 'COMPLETED', 'item': item, 'pairs': len(features),
        'output_sha256': retained, 'submission_ready': False, 'submission_count': 0})
    print(json.dumps({'index': args.index, 'state': 'COMPLETED', 'pairs': len(features)}), flush=True)


def load_split(args, protocol, split):
    rows, baseline_cohorts, hashes = [], [], {}
    for item in protocol['cohorts']:
        if item['stage'] != split:
            continue
        directory = args.out_root / f"cohort_{item['index']:03d}"
        result = json.loads((directory / 'result.json').read_text())
        for name in ['features.json', 'truth.json', 'baseline_predictions.json', 'identifiers.json']:
            path = (directory / name).resolve()
            expected = result['output_sha256'][str(path)]
            if sha256(path) != expected:
                raise ValueError('Generated panel output changed')
            hashes[str(path)] = expected
        truth = json.loads((directory / 'truth.json').read_text())
        target = {(row['indID1'], row['indID2']): row['relation'] for row in truth}
        features = json.loads((directory / 'features.json').read_text())
        ids = sorted(json.loads((directory / 'identifiers.json').read_text()).values())
        if set(target) != set(itertools.combinations(ids, 2)):
            raise ValueError('Truth omits roster pairs')
        if {(row['indID1'], row['indID2']) for row in features} != set(target) or len(features) != len(target):
            raise ValueError('Incomplete detected feature coverage')
        cohort_id = f"detected_panel_{item['index']:03d}"
        rows.extend({**row, 'relation': target[(row['indID1'], row['indID2'])], 'cohort_id': cohort_id} for row in features)
        baseline_cohorts.append({'cohort_id': cohort_id, 'truth': truth,
            'predictions': json.loads((directory / 'baseline_predictions.json').read_text())})
    return rows, baseline_cohorts, hashes


def fit(args, protocol):
    from sklearn.ensemble import ExtraTreesClassifier
    output = args.out_root / 'model_evaluation'
    output.mkdir(parents=True, exist_ok=False)
    hashes = bound_inputs(args)
    training, _, bound = load_split(args, protocol, 'training')
    hashes.update(bound)
    names = sorted(set(training[0]) - {'indID1', 'indID2', 'relation', 'cohort_id'})
    bundles = {}
    for variant, feature_names in [('context', names), ('pair_only', [name for name in names if not name.startswith('context_')])]:
        model = ExtraTreesClassifier(n_estimators=500, min_samples_leaf=5, class_weight='balanced',
            max_features=1., random_state=protocol['model_seed'], n_jobs=4)
        model.fit([[row[name] for name in feature_names] for row in training], [row['relation'] for row in training])
        bundle = {'model': model, 'feature_names': feature_names, 'candidate': protocol['candidate'],
            'variant': variant, 'protocol_sha256': sha256(args.protocol), 'submission_ready': False}
        bundles[variant] = bundle
        joblib.dump(bundle, output / f'{variant}_model.joblib')
    calibration, baseline, bound = load_split(args, protocol, 'calibration')
    hashes.update(bound)
    selected, scores = {}, {}
    labels = list(bundles['context']['model'].classes_)
    target = np.asarray([labels.index(row['relation']) for row in calibration])
    counts = np.bincount(target, minlength=len(labels))
    if np.any(counts == 0):
        raise ValueError('Calibration lacks a task label')
    weights = 1 / counts[target]
    for variant, bundle in bundles.items():
        model = bundle['model']
        values = model.predict_proba([[row[name] for name in bundle['feature_names']] for row in calibration])
        losses = []
        for value in protocol['scalar_temperature_grid']:
            adjusted = temperature(values, value)
            loss = -float(np.sum(weights * np.log(np.maximum(adjusted[np.arange(len(values)), target], 1e-12))) / weights.sum())
            losses.append({'temperature': value, 'balanced_log_loss': loss})
        selected[variant] = min(losses, key=lambda row: (row['balanced_log_loss'], row['temperature']))
        scorer_input, score = score_family_rows(calibration, temperature(values, selected[variant]['temperature']), model.classes_)
        save_json(output / f'{variant}_calibration_scorer_input.json', scorer_input)
        scores[f'{variant}_calibration'] = score
    freeze = {'hashes': hashes, 'temperatures': selected, 'model_hashes': {variant: sha256(output / f'{variant}_model.joblib') for variant in bundles},
        'candidate': protocol['candidate'], 'development_outcomes_not_yet_read': True}
    save_json(output / 'model_calibration_freeze.json', freeze)
    development, baseline, bound = load_split(args, protocol, 'development')
    hashes.update(bound)
    for variant, bundle in bundles.items():
        model = bundle['model']
        values = model.predict_proba([[row[name] for name in bundle['feature_names']] for row in development])
        values = temperature(values, selected[variant]['temperature'])
        scorer_input, score = score_family_rows(development, values, model.classes_)
        save_json(output / f'{variant}_development_scorer_input.json', scorer_input)
        save_json(output / f'{variant}_development_predictions.json', [{**row, 'cohort_id': cohort['cohort_id']}
            for cohort in scorer_input['cohorts'] for row in cohort['predictions']])
        scores[f'{variant}_development'] = metrics(scorer_input['cohorts'], labels)
    baseline_score = metrics(baseline, labels)
    save_json(output / 'baseline_development_scorer_input.json', {'labels': labels, 'cohorts': baseline})
    score = scores['context_development']
    gates = {'all_labels': score['complete_label_coverage'], 'exact_accuracy': score['exact_accuracy'] >= .8,
        'macro_recall': score['macro_recall'] >= .7, 'minimum_label_recall': score['minimum_label_recall'] >= .5,
        'gain_over_baseline': score['exact_accuracy'] - baseline_score['exact_accuracy'] >= .1,
        'no_cohort_regression': all(value >= baseline_score['cohort_accuracy'][cohort] for cohort, value in score['cohort_accuracy'].items())}
    verify_inputs(hashes)
    save_json(output / 'result.json', {'stage': 'DEVELOPMENT', 'scores': scores, 'baseline': baseline_score,
        'gates': gates, 'advance_to_confirmation': all(gates.values()), 'input_sha256': hashes,
        'submission_ready': False, 'submission_count': 0})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['generate', 'fit'])
    parser.add_argument('--asset-root', type=Path, required=True)
    parser.add_argument('--source-root', type=Path, required=True)
    parser.add_argument('--phase-model', type=Path, required=True)
    parser.add_argument('--protocol', type=Path, required=True)
    parser.add_argument('--out-root', type=Path, required=True)
    parser.add_argument('--index', type=int)
    args = parser.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):
        parser.error('Use Slurm for the detected training panel')
    verify_inputs(json.loads((args.source_root / 'source_manifest.json').read_text()))
    protocol = json.loads(args.protocol.read_text())
    if protocol['candidate'] != 'detected-multiway-family-forest-v2':
        raise ValueError('Unexpected policy')
    if sha256(args.phase_model) != protocol['frozen_sibling_model_sha256']:
        raise ValueError('Frozen sibling detector changed')
    if args.command == 'generate':
        generate(args, protocol)
    else:
        fit(args, protocol)


if __name__ == '__main__':
    main()
