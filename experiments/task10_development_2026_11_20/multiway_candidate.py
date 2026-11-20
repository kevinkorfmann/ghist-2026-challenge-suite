"""Prospectively fixed development test of inferred sibling pools and multiway IBD."""
from __future__ import annotations

import argparse
import itertools
import json
import os
import subprocess
from collections import defaultdict
from pathlib import Path

import numpy as np

from ghist_relatedness.labels import HUMAN_LABELS
from scripts.develop_task10_detected_ibd import (
    family_truth,
    file_hashes,
    genome_length,
    merge_spans,
    parents_from_fam,
    phase_pair_features,
    read_exact_phase,
    read_map,
    save_json,
    score_predictions,
    sha256,
    verify_inputs,
)


def topology_design(count, seed, excluded=()):
    rng = np.random.default_rng(seed)
    used, result = set(map(tuple, excluded)), []
    while len(result) < count:
        lower = 1 if len(result) % 2 == 0 else 2
        depths = tuple(sorted(map(int, rng.integers(lower, lower + 4, size=int(rng.integers(4, 9))))))
        if depths not in used:
            used.add(depths)
            result.append(depths)
    return result


def whole_family_definition(topologies, prefix):
    lines, families, expected = [], {}, {}
    for index, depths in enumerate(topologies):
        name = f'{prefix}{index:04d}'
        lines.append(f'def {name} 1 {max(depths) + 1}')
        for generation in range(2, max(depths) + 2):
            omitted = [f'{branch}n' for branch, depth in enumerate(depths, 1) if depth != generation - 1]
            lines.append(f'{generation} 2 {len(depths)} ' + ' '.join(omitted))
        samples = {f'{name}1_g{depth + 1}-b{branch}-i{child}': (branch, depth)
            for branch, depth in enumerate(depths, 1) for child in [1, 2]}
        truth = {}
        for a, b in itertools.combinations(sorted(samples), 2):
            branch_a, depth_a = samples[a]
            branch_b, depth_b = samples[b]
            relation = (1, 1, 2) if branch_a == branch_b else (*sorted([depth_a, depth_b]), 2)
            truth[(a, b)] = HUMAN_LABELS[relation]
        families[name] = truth
        expected.update(truth)
    return '\n'.join(lines) + '\n', families, expected


def inferred_groups(samples, features, phase_model):
    model = phase_model['model']
    probabilities = model.predict_proba([[row[key] for key in phase_model['feature_names']] for row in features])
    sibling_column = list(model.classes_).index('siblings')
    adjacency = {sample: set() for sample in samples}
    for row, values in zip(features, probabilities, strict=True):
        if values[sibling_column] >= .95 and row['ibd2_fraction_5'] >= .10:
            a, b = row['indID1'], row['indID2']
            adjacency[a].add(b)
            adjacency[b].add(a)
    remaining, groups = set(samples), []
    while remaining:
        pending, component = [min(remaining)], set()
        while pending:
            sample = pending.pop()
            if sample not in component:
                component.add(sample)
                pending.extend(adjacency[sample] - component)
        remaining -= component
        if all(component - {sample} <= adjacency[sample] for sample in component):
            groups.append(tuple(sorted(component)))
        else:
            groups.extend((sample,) for sample in sorted(component))
    return sorted(groups), probabilities


def component_integrals(groups, segments, genetic_map, genome_cm, threshold):
    """Count distinct local shared homolog components, excluding impossible merges."""
    membership = {sample: index for index, group in enumerate(groups) for sample in group}
    events = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    for a, b in itertools.combinations(sorted(membership), 2):
        for chrom, rows in segments.get((a, b), {}).items():
            by_edge = defaultdict(list)
            for start, end, h1, h2 in rows:
                by_edge[((a, h1), (b, h2))].append((start, end))
            for edge, spans in by_edge.items():
                for start, end in merge_spans(spans):
                    lo, hi = map(float, np.interp([start, end], *genetic_map[chrom]))
                    if hi - lo >= threshold:
                        events[chrom][lo][edge] += 1
                        events[chrom][hi][edge] -= 1
    totals = dict.fromkeys(['union', 'copies', 'multiple', 'conflict', 'complete'], 0.)
    required_types = set(itertools.combinations(range(len(groups)), 2))
    for chromosome_events in events.values():
        active, previous = defaultdict(int), None
        for position in sorted(chromosome_events):
            edges = [edge for edge, count in active.items() if count > 0]
            if previous is not None and edges:
                parent = {node: node for edge in edges for node in edge}

                def find(node):
                    while parent[node] != node:
                        parent[node] = parent[parent[node]]
                        node = parent[node]
                    return node

                for a, b in edges:
                    parent[find(a)] = find(b)
                components, edge_types = defaultdict(set), defaultdict(set)
                for node in parent:
                    components[find(node)].add(node)
                for a, b in edges:
                    if membership[a[0]] != membership[b[0]]:
                        edge_types[find(a)].add(tuple(sorted([membership[a[0]], membership[b[0]]])))
                shared, complete, conflicts = 0, 0, False
                for root, nodes in components.items():
                    people = [sample for sample, _ in nodes]
                    if len(people) != len(set(people)):
                        conflicts = True
                        continue
                    if len({membership[sample] for sample in people}) == len(groups):
                        shared += 1
                        complete += required_types <= edge_types[root]
                length = position - previous
                totals['union'] += length * (shared > 0)
                totals['copies'] += length * shared
                totals['multiple'] += length * (shared >= 2)
                totals['conflict'] += length * conflicts
                totals['complete'] += length * (complete > 0)
            for edge, count in chromosome_events[position].items():
                active[edge] += count
            previous = position
    return {name: value / genome_cm for name, value in totals.items()}


def family_features(samples, segments, genetic_map, genome_cm, phase_model):
    rows = [{'indID1': a, 'indID2': b, **phase_pair_features(segments.get((a, b), {}), genetic_map, genome_cm)}
        for a, b in itertools.combinations(sorted(samples), 2)]
    groups, frozen_probabilities = inferred_groups(samples, rows, phase_model)
    membership = {sample: index for index, group in enumerate(groups) for sample in group}
    cached = {}

    def integral(indices, threshold):
        key = (tuple(sorted(indices)), threshold)
        if key not in cached:
            cached[key] = component_integrals([groups[index] for index in key[0]], segments,
                genetic_map, genome_cm, threshold)
        return cached[key]

    for row in rows:
        a, b = membership[row['indID1']], membership[row['indID2']]
        low, high = sorted([len(groups[a]), len(groups[b])])
        row.update(context_same_group=float(a == b), context_group_size_low=float(low),
            context_group_size_high=float(high), context_groups=float(len(groups)))
        for threshold in [5, 10]:
            pooled = integral([a, b], threshold) if a != b else dict.fromkeys(['union', 'copies', 'multiple', 'conflict', 'complete'], 0.)
            for name, value in pooled.items():
                row[f'context_pool_{name}_{threshold}'] = value
            triples = [integral([a, b, c], threshold) for c in range(len(groups)) if c not in [a, b]] if a != b else []
            for name in ['union', 'copies', 'multiple', 'conflict', 'complete', 'relative_union']:
                values = [triple[name] if name != 'relative_union' else
                    triple['union'] / max(pooled['union'], 1e-12) for triple in triples]
                statistics = [float(np.mean(values)), *map(float, np.quantile(values, [0., .25, .5, .75, 1.]))] if values else [0.] * 6
                for statistic, value in zip(['mean', 'min', 'q25', 'median', 'q75', 'max'], statistics, strict=True):
                    row[f'context_triple_{name}_{statistic}_{threshold}'] = value
    return rows, groups, frozen_probabilities


def score_family_rows(rows, probabilities, classes):
    cohorts = {}
    for row, values in zip(rows, probabilities, strict=True):
        cohort = cohorts.setdefault(row['cohort_id'], {'cohort_id': row['cohort_id'], 'truth': [], 'predictions': []})
        pair = {key: row[key] for key in ['indID1', 'indID2']}
        cohort['truth'].append({**pair, 'relation': row['relation']})
        cohort['predictions'].append({**pair, 'relation': str(classes[int(np.argmax(values))]),
            'probabilities': dict(zip(map(str, classes), map(float, values), strict=True))})
    def distinct_rows(kind):
        return [{**row, 'indID1': f"{cohort['cohort_id']}::{row['indID1']}",
            'indID2': f"{cohort['cohort_id']}::{row['indID2']}"}
            for cohort in cohorts.values() for row in cohort[kind]]
    truth, predictions = distinct_rows('truth'), distinct_rows('predictions')
    score = score_predictions(truth, predictions)
    recalls = [entry['correct'] / entry['pairs'] for entry in score['strata'].values()]
    score['macro_recall'] = float(np.mean(recalls))
    score['minimum_label_recall'] = min(recalls)
    score['complete_label_coverage'] = set(score['strata']) == set(HUMAN_LABELS.values())
    return {'labels': list(HUMAN_LABELS.values()), 'cohorts': list(cohorts.values())}, score


def main():
    import joblib
    from sklearn.ensemble import ExtraTreesClassifier

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--asset-root', type=Path, required=True)
    parser.add_argument('--source-root', type=Path, required=True)
    parser.add_argument('--phase-model', type=Path, required=True)
    parser.add_argument('--protocol', type=Path, required=True)
    parser.add_argument('--out-dir', type=Path, required=True)
    args = parser.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):
        parser.error('Full simulation and model fitting require Slurm')
    assets, output = args.asset_root.resolve(), args.out_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    protocol = json.loads(args.protocol.read_text())
    if protocol['candidate'] != 'multiway-family-forest-v1' or sha256(args.phase_model) != protocol['frozen_pair_model_sha256']:
        raise ValueError('Unexpected frozen candidate or model')
    phase_model = joblib.load(args.phase_model)
    pedsim = assets / 'external/ped-sim/ped-sim'
    simmap = assets / 'external/ped-sim-map/refined_mf.simmap'
    map_path = assets / 'external/ped-sim-map/refined_average.plink.map'
    interference = assets / 'external/ped-sim/interfere/nu_p_campbell.tsv'
    hashes = file_hashes(args.source_root)
    hashes.update({str(path.resolve()): sha256(path) for path in [Path(__file__), args.phase_model,
        args.protocol, pedsim, simmap, map_path, interference]})
    design = {}
    for stage in ['training', 'exact_diagnostic']:
        spec = protocol[stage]
        design[stage] = topology_design(spec['families'], spec['topology_seed'], design.get('training', []))
    save_json(output / 'preregistered_design.json', {'protocol': protocol, 'topologies': design, 'hashes': hashes})
    genetic_map, genome_cm = read_map(map_path), genome_length(map_path)
    datasets, frozen_outputs, group_diagnostics = {}, {}, {}
    dropout_rng = np.random.default_rng(protocol['dropout_seed'])
    for stage in ['training', 'exact_diagnostic']:
        definition, families, expected = whole_family_definition(design[stage], stage)
        definition_path = output / f'{stage}.def'
        definition_path.write_text(definition)
        with (output / f'{stage}.log').open('w') as log:
            subprocess.run([str(pedsim), '-d', str(definition_path), '-m', str(simmap),
                '-o', str(output / stage), '--intf', str(interference), '--seed', str(protocol[stage]['meiosis_seed']),
                '--seg_haps', '--fam'], check=True, stdout=log, stderr=subprocess.STDOUT)
        parents = parents_from_fam(output / f'{stage}-everyone.fam')
        for truth in families.values():
            samples = {sample for pair in truth for sample in pair}
            observed = {(row['indID1'], row['indID2']): row['relation']
                for row in family_truth(parents, {sample: sample for sample in samples})}
            if observed != truth:
                raise ValueError('Simulated pedigree parentage differs from declared topology')
        segments = read_exact_phase(output / f'{stage}.seg', expected)
        stages = [stage, 'dropout_diagnostic'] if stage == 'exact_diagnostic' else [stage]
        for feature_stage in stages:
            all_rows, all_frozen, diagnostics = [], [], []
            for family, truth in families.items():
                samples = sorted({sample for pair in truth for sample in pair})
                subset = {pair: segments[pair] for pair in truth if pair in segments}
                if feature_stage == 'dropout_diagnostic':
                    subset = {pair: {chrom: [row for row in rows if dropout_rng.random() >= .2]
                        for chrom, rows in chromosomes.items()} for pair, chromosomes in subset.items()}
                rows, groups, probabilities = family_features(samples, subset, genetic_map, genome_cm, phase_model)
                for row in rows:
                    row['relation'] = truth[(row['indID1'], row['indID2'])]
                    row['cohort_id'] = family
                all_rows.extend(rows)
                all_frozen.extend(probabilities)
                diagnostics.append({'cohort_id': family, 'groups': groups})
            datasets[feature_stage], frozen_outputs[feature_stage] = all_rows, np.asarray(all_frozen)
            group_diagnostics[feature_stage] = diagnostics
            save_json(output / f'{feature_stage}_features.json', all_rows)
            print(feature_stage, len(all_rows), flush=True)
    save_json(output / 'inferred_groups.json', group_diagnostics)
    names = sorted(set(datasets['training'][0]) - {'indID1', 'indID2', 'cohort_id', 'relation'})
    summary = {}
    for variant, feature_names in [('context', names), ('pair_only', [name for name in names if not name.startswith('context_')])]:
        classifier = ExtraTreesClassifier(n_estimators=500, min_samples_leaf=5, class_weight='balanced',
            random_state=protocol['model_seed'], max_features=1., n_jobs=4)
        classifier.fit([[row[name] for name in feature_names] for row in datasets['training']],
            [row['relation'] for row in datasets['training']])
        joblib.dump({'model': classifier, 'feature_names': feature_names, 'protocol_sha256': sha256(args.protocol),
            'candidate': protocol['candidate'], 'variant': variant, 'submission_ready': False}, output / f'{variant}_model.joblib')
        for stage in ['exact_diagnostic', 'dropout_diagnostic']:
            probabilities = classifier.predict_proba([[row[name] for name in feature_names] for row in datasets[stage]])
            scorer_input, score = score_family_rows(datasets[stage], probabilities, classifier.classes_)
            save_json(output / f'{variant}_{stage}_scorer_input.json', scorer_input)
            summary[f'{variant}_{stage}'] = score
    for stage in ['exact_diagnostic', 'dropout_diagnostic']:
        scorer_input, score = score_family_rows(datasets[stage], frozen_outputs[stage], phase_model['model'].classes_)
        save_json(output / f'frozen_{stage}_scorer_input.json', scorer_input)
        summary[f'frozen_{stage}'] = score
    context, pair = summary['context_exact_diagnostic'], summary['pair_only_exact_diagnostic']
    screen = context['exact_accuracy'] - pair['exact_accuracy'] >= .05 and context['macro_recall'] >= pair['macro_recall']
    verify_inputs(hashes)
    save_json(output / 'result.json', {'stage': 'DEVELOPMENT', 'scores': summary, 'exact_screen_passed': screen,
        'model_hashes': {variant: sha256(output / f'{variant}_model.joblib') for variant in ['context', 'pair_only']},
        'submission_ready': False, 'submission_count': 0})


if __name__ == '__main__':
    main()
