"""Separate Task 10 development experiment; never promotes a production release."""
from __future__ import annotations

import argparse
import csv
import gzip
import itertools
import json
import os
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

try:
    from scripts.fit_human_family_relationships import summarize_raw_lengths
    from scripts.reconstruct_bonsai_pedigree import (
        coverage_intervals,
        read_haplotype_segments,
        read_map,
    )
    from scripts.summarize_hapibd import genome_length
    from scripts.validate_task10_complete_pipeline import (
        family_truth,
        frozen_pair_predictions,
        map_bounds,
        parents_from_fam,
        save_json,
        score_predictions,
        sha256,
    )
except ModuleNotFoundError:
    from fit_human_family_relationships import summarize_raw_lengths
    from reconstruct_bonsai_pedigree import (
        coverage_intervals,
        read_haplotype_segments,
        read_map,
    )
    from summarize_hapibd import genome_length
    from validate_task10_complete_pipeline import (
        family_truth,
        frozen_pair_predictions,
        map_bounds,
        parents_from_fam,
        save_json,
        score_predictions,
        sha256,
    )


def complete_features(segments: dict, samples: list[str], genome_cm: float) -> list[dict]:
    """All roster pairs including zero detected sharing, independent of scoring truth."""
    if len(set(samples)) != len(samples) or len(samples) < 2:
        raise ValueError("The genotype roster must contain unique sample identifiers")
    valid = set(itertools.combinations(sorted(samples), 2))
    if set(segments) - valid:
        raise ValueError("Segments reference a pair outside the genotype roster")
    return [dict(indID1=a, indID2=b, **summarize_raw_lengths(segments.get((a, b), []), genome_cm))
            for a, b in sorted(valid)]


def genotype_roster(path: Path) -> list[str]:
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as handle:
        for line in handle:
            if line.startswith('#CHROM'):
                return line.rstrip().split('\t')[9:]
    raise ValueError('VCF has no sample header')


def detected_segments(path: Path, map_path: Path, representation: str) -> dict:
    grouped, _ = read_haplotype_segments(path)
    genetic_map = read_map(map_path)
    result = defaultdict(list)
    for (a, b, chrom), intervals in grouped.items():
        bp, cm = genetic_map[chrom]
        if representation == 'diploid_state':
            pieces = coverage_intervals(intervals)
        elif representation == 'raw_haplotype':
            pieces = [(start, end, False) for start, end, _, _ in intervals]
        else:
            raise ValueError('Unknown feature representation')
        for start, end, ibd2 in pieces:
            if start < bp[0] or end > bp[-1]:
                raise ValueError('Segment exceeds map bounds')
            values = np.interp([start, end], bp, cm)
            length = float(values[1] - values[0])
            if length > 0:
                result[(a, b)].append((chrom, length))
                if ibd2:
                    result[(a, b)].append((chrom, length))
    return result


def exact_segments(path: Path, identifiers: dict) -> dict:
    result = defaultdict(list)
    for line in path.read_text().splitlines():
        a, b, chrom, _, _, state, _, _, length = line.split()
        if a not in identifiers or b not in identifiers:
            continue
        pair = tuple(sorted((identifiers[a], identifiers[b])))
        result[pair].append((chrom, float(length)))
        if state == 'IBD2':
            result[pair].append((chrom, float(length)))
    return result


def write_features(rows: list[dict], path: Path) -> None:
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def evaluate_representations(cohort: Path, assets: Path, output: Path) -> dict:
    output.mkdir(parents=True, exist_ok=False)
    map_path = assets / 'external/ped-sim-map/refined_average.plink.map'
    identifiers = json.loads((cohort / 'identifiers.json').read_text())
    roster = genotype_roster(cohort / 'family.normalized.vcf.gz')
    truth = json.loads((cohort / 'truth.json').read_text())
    sources = {
        'exact_diploid_oracle': exact_segments(cohort / 'family.seg', identifiers),
        'raw_haplotype': detected_segments(cohort / 'detected.ibd.gz', map_path, 'raw_haplotype'),
        'diploid_state': detected_segments(cohort / 'detected.ibd.gz', map_path, 'diploid_state'),
    }
    candidate = assets / 'artifacts/task10_large_family/clean_spectrum_v1/model.joblib'
    baseline = assets / 'artifacts/task10_large_family/pair_model_min5/model.joblib'
    results = {}
    for name, segments in sources.items():
        rows = complete_features(segments, roster, genome_length(map_path))
        feature_path = output / f'{name}.csv'
        write_features(rows, feature_path)
        predictions = frozen_pair_predictions(feature_path, candidate, baseline, output / f'{name}.json')
        results[name] = score_predictions(truth, predictions)
    save_json(output / 'result.json', {'development_only': True, 'results': results,
        'cohort': str(cohort), 'input_hashes': {str(cohort / name): sha256(cohort / name)
        for name in ['identifiers.json', 'truth.json', 'family.seg', 'detected.ibd.gz']},
        'submission_ready': False, 'submission_count': 0})
    return results



TOPOLOGIES = {
    'cousins': {'depths': [2, 2, 3, 3, 4, 4, 5, 5], 'samples_per_branch': 2},
    'avuncular': {'depths': [1, 2, 3, 4], 'samples_per_branch': 2},
}


def pedigree_definition(topology: str) -> str:
    depths = TOPOLOGIES[topology]['depths']
    branches = len(depths)
    lines = [f'def validation 1 {max(depths) + 1}']
    for generation in range(2, max(depths) + 2):
        omitted = [str(i + 1) + 'n' for i, depth in enumerate(depths) if depth != generation - 1]
        lines.append(f"{generation} 2 {branches} " + ' '.join(omitted))
    return '\n'.join(lines) + '\n'



def derived_seed(cohort_seed: int, chromosome: int, stage: int) -> int:
    """Disjoint identified streams; stage1 ancestry,2 mutation,3 meioses,4 dropout."""
    value = int(np.random.SeedSequence([cohort_seed, chromosome, stage]).generate_state(1, dtype=np.uint32)[0])
    return value % (2**32 - 1) + 1

def simulated_founder_chromosome(job: tuple) -> dict:
    import msprime
    chrom, start, end, cm, seeds, count, filename = job
    seed, mutation_seed = seeds
    length = end - start + 1
    ancestry = msprime.sim_ancestry(samples=count, population_size=10000,
        sequence_length=length, recombination_rate=cm / (100 * length), random_seed=seed)
    mutated = msprime.sim_mutations(ancestry, rate=1.25e-8,
        model=msprime.BinaryMutationModel(), random_seed=mutation_seed)
    retained, last = 0, -1000
    with gzip.open(filename, 'wt', compresslevel=1) as handle:
        for variant in mutated.variants():
            position = int(variant.site.position) + start
            minor = min(int(variant.genotypes.sum()), 2 * count - int(variant.genotypes.sum()))
            if position - last < 500 or minor / (2 * count) < 2 / 28:
                continue
            calls = '\t'.join(f'{a}|{b}' for a, b in np.asarray(variant.genotypes).reshape(count, 2))
            handle.write(f'{chrom}\t{position}\t.\tA\tC\t.\tPASS\t.\tGT\t{calls}\n')
            retained, last = retained + 1, position
    return {'chromosome': chrom, 'markers': retained, 'ancestry_seed': seed,
            'mutation_seed': mutation_seed}


def generate_founders(bounds: dict, seed: int, count: int, root: Path, workers: int) -> Path:
    root.mkdir()
    jobs = [(chrom, *bounds[chrom], (derived_seed(seed, int(chrom), 1), derived_seed(seed, int(chrom), 2)), count,
             str(root / f'chr{chrom}.vcf.gz')) for chrom in sorted(bounds, key=int)]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        records = list(pool.map(simulated_founder_chromosome, jobs))
    output = root / 'founders.vcf.gz'
    with gzip.open(output, 'wt', compresslevel=1) as writer:
        writer.write('##fileformat=VCFv4.2\n')
        writer.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        for chrom in sorted(bounds, key=int):
            writer.write(f'##contig=<ID={chrom},length={bounds[chrom][1] + 1}>\n')
        writer.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t')
        writer.write('\t'.join(f'founder_{i}' for i in range(count)) + '\n')
        for job in jobs:
            with gzip.open(job[-1], 'rt') as reader:
                shutil.copyfileobj(reader, writer)
    save_json(root / 'generation.json', {'founders': count, 'chromosomes': records})
    return output


def select_vcf(source: Path, destination: Path, topology: str, parents: dict) -> dict:
    depths = TOPOLOGIES[topology]['depths']
    identifiers = {}
    indices = []
    variants = 0
    with gzip.open(source, 'rt') as reader, gzip.open(destination, 'wt', compresslevel=1) as writer:
        for line in reader:
            if line.startswith('#CHROM'):
                fields = line.rstrip().split('\t')
                for i, sample in enumerate(fields[9:]):
                    match = re.fullmatch(r'validation1_g(\d+)-b(\d+)-i([12])', sample)
                    if match and len(parents[sample]) == 2:
                        generation, branch, _ = map(int, match.groups())
                        if generation - 1 == depths[branch - 1]:
                            identifiers[sample] = f'sample_{len(identifiers)}'
                            indices.append(i + 9)
                writer.write('\t'.join(fields[:9] + list(identifiers.values())) + '\n')
            elif line.startswith('#'):
                writer.write(line)
            else:
                fields = line.rstrip().split('\t')
                calls = [fields[i] for i in indices]
                if any(call.split(':')[0] not in {'0|0', '0|1', '1|0', '1|1'} for call in calls):
                    raise ValueError('Phased, complete synthetic genotype contract failed')
                writer.write('\t'.join(fields[:9] + calls) + '\n')
                variants += 1
    if len(identifiers) != 2 * len(depths) or not variants:
        raise ValueError(f'Unexpected genotype roster: {identifiers}')
    return identifiers


def file_hashes(root: Path) -> dict:
    return {str(path.resolve()): sha256(path) for folder in ['scripts', 'src']
            for path in sorted((root / folder).rglob('*.py'))}


def verify_inputs(hashes: dict) -> None:
    for path, expected in hashes.items():
        if sha256(Path(path)) != expected:
            raise RuntimeError(f'Immutable input changed: {path}')


def run_cohort(args) -> None:
    root, assets, output = args.project_root.resolve(), args.asset_root.resolve(), args.out_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    protocol = json.loads(args.protocol.read_text())
    item = next((row for row in protocol['cohorts'] if row['index'] == args.index), None)
    if item is None:
        raise ValueError('Cohort must be in the preregistered protocol')
    if protocol['candidate'] != 'frozen-classifier-diploid-state-whole-family-v1':
        raise ValueError('Unrecognized candidate')
    topology, seed = item['topology'], item['seed']
    env = {**os.environ, 'PYTHONPATH': f'{root / "src"}:{root / "scripts"}',
           'PYTHONDONTWRITEBYTECODE': '1'}
    commands = []
    def command(name, words):
        with (output / f'{name}.log').open('w') as log:
            result = subprocess.run(list(map(str, words)), check=False, cwd=root, env=env,
                stdout=log, stderr=subprocess.STDOUT)
        commands.append({'stage': name, 'command': list(map(str, words)), 'exit_code': result.returncode})
        save_json(output / 'commands.json', commands)
        print(name, result.returncode, flush=True)
        if result.returncode:
            raise RuntimeError(f'{name} failed')
    def script(name):
        return [sys.executable, root / 'scripts' / name]
    map_path = assets / 'external/ped-sim-map/refined_average.plink.map'
    simmap = assets / 'external/ped-sim-map/refined_mf.simmap'
    pedsim = assets / 'external/ped-sim/ped-sim'
    interference = assets / 'external/ped-sim/interfere/nu_p_campbell.tsv'
    detector = assets / 'external/hap-ibd/hap-ibd.jar'
    candidate = assets / 'artifacts/task10_large_family/clean_spectrum_v1/model.joblib'
    baseline = assets / 'artifacts/task10_large_family/pair_model_min5/model.joblib'
    hashes = file_hashes(root)
    hashes.update({str(path): sha256(path) for path in
        [map_path, simmap, pedsim, interference, detector, candidate, baseline, args.protocol.resolve()]})
    save_json(output / 'preregistered_design.json', {'item': item, 'protocol': protocol,
        'hashes': hashes, 'submission_count': 0, 'scientific_release_qualified': False,
        'derived_seeds': {'chromosomes': {chrom: {'ancestry': derived_seed(seed, int(chrom), 1),
            'mutation': derived_seed(seed, int(chrom), 2), 'dropout_unused': derived_seed(seed, int(chrom), 4)}
            for chrom in map_bounds(map_path)}, 'meioses': derived_seed(seed, 0, 3)}})
    depths = TOPOLOGIES[topology]['depths']
    founder_count = 2 + len(depths) * (max(depths) - 1)
    founders = generate_founders(map_bounds(map_path), seed, founder_count, output / 'founders', args.workers)
    definition = output / 'family.def'
    definition.write_text(pedigree_definition(topology))
    command('pedsim', [pedsim, '-d', definition, '-m', simmap, '-i', founders,
        '-o', output / 'family', '--intf', interference, '--seed', derived_seed(seed, 0, 3),
        '--fam', '--keep_phase', '--err_rate', 0, '--miss_rate', 0])
    parents = parents_from_fam(output / 'family-everyone.fam')
    normalized = output / 'family.normalized.vcf.gz'
    identifiers = select_vcf(output / 'family.vcf.gz', normalized, topology, parents)
    save_json(output / 'identifiers.json', identifiers)
    truth = family_truth(parents, identifiers)
    save_json(output / 'truth.json', truth)
    hashes.update({str(output / name): sha256(output / name) for name in
                   ['truth.json', 'identifiers.json', 'family.normalized.vcf.gz', 'family-everyone.fam']})
    save_json(output / 'inference_input_hashes.json', hashes)
    prepared = output / 'prepared.vcf.gz'
    command('prepare', script('prepare_hapibd_vcf.py') + ['--input', normalized, '--output', prepared])
    command('detector', ['java', '-Xmx14g', '-jar', detector, f'gt={prepared}', f'map={map_path}',
        f'out={output / "detected"}', 'min-seed=1.0', 'min-extend=0.2', 'min-output=2.0',
        'min-markers=100', 'max-gap=1000', f'nthreads={args.workers}'])
    results = evaluate_representations(output, assets, output / 'representation_comparison')
    clusters = output / 'clusters.json'
    save_json(clusters, {'families': [sorted(identifiers.values())], 'samples': len(identifiers),
                         'family_sizes': [len(identifiers)], 'source': 'Task10 single-family contract'})
    for name, degree, types in [('fast', 5, 2), ('medium', 5, 3)]:
        command(name, script('reconstruct_human_branch_pedigrees.py') + ['--ibd', output / 'detected.ibd.gz',
            '--map', map_path, '--clusters', clusters, '--threshold-cm', 10,
            '--max-add-degree', degree, '--max-rel-types', types, '--out', output / f'{name}.json'])
    command('consensus', script('consensus_human_family_pedigree.py') + ['--pair-predictions',
        output / 'representation_comparison/diploid_state.json', '--bonsai', output / 'fast.json',
        output / 'medium.json', '--topology-boost', 1, '--require-flat-degree-agreement',
        '--out-dir', output / 'consensus'])
    with (output / 'consensus/candidate.csv').open() as handle:
        predictions = list(csv.DictReader(handle))
    scores = score_predictions(truth, predictions)
    verify_inputs(hashes)
    cohort = {'cohort_id': f'{topology}_{seed}', 'truth': truth, 'predictions': predictions}
    save_json(output / 'scoring_cohort.json', cohort)
    retained = {str(path.relative_to(output)): sha256(path) for path in output.rglob('*')
        if path.is_file() and path.name not in {'result.json'}}
    save_json(output / 'result.json', {'state': 'COMPLETED', 'item': item, 'candidate': scores,
        'comparisons': results, 'outputs': retained,
        'protocol_sha256': sha256(args.protocol), 'submission_ready': False, 'submission_count': 0})
    print(json.dumps({'candidate': scores, 'comparisons': results}), flush=True)



def run_whole_family(args) -> None:
    root, assets, cohort, output = args.project_root.resolve(), args.asset_root.resolve(), args.cohort.resolve(), args.out_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    roster = genotype_roster(cohort / 'family.normalized.vcf.gz')
    map_path = assets / 'external/ped-sim-map/refined_average.plink.map'
    truth = json.loads((cohort / 'truth.json').read_text())
    features = output / 'features.csv'
    write_features(complete_features(detected_segments(cohort / 'detected.ibd.gz', map_path,
        'diploid_state'), roster, genome_length(map_path)), features)
    predicted = frozen_pair_predictions(features,
        assets / 'artifacts/task10_large_family/clean_spectrum_v1/model.joblib',
        assets / 'artifacts/task10_large_family/pair_model_min5/model.joblib', output / 'pair_predictions.json')
    clusters = output / 'clusters.json'
    save_json(clusters, {'families': [sorted(roster)], 'samples': len(roster),
        'family_sizes': [len(roster)], 'source': 'Task10 single-family contract'})
    env = {**os.environ, 'PYTHONPATH': f'{root / "src"}:{root / "scripts"}', 'PYTHONDONTWRITEBYTECODE': '1'}
    hashes = file_hashes(root)
    hashes.update({str(cohort / name): sha256(cohort / name) for name in
        ['truth.json', 'detected.ibd.gz', 'family.normalized.vcf.gz']})
    save_json(output / 'protocol.json', {'design': 'whole-family development diagnosis',
        'settings': [[5, 2], [5, 3]], 'hashes': hashes, 'submission_ready': False})
    for name, degree, types in [('fast', 5, 2), ('medium', 5, 3)]:
        with (output / f'{name}.log').open('w') as log:
            subprocess.run([sys.executable, root / 'scripts/reconstruct_human_branch_pedigrees.py',
                '--ibd', cohort / 'detected.ibd.gz', '--map', map_path, '--clusters', clusters,
                '--threshold-cm', '10', '--max-add-degree', str(degree), '--max-rel-types', str(types),
                '--out', output / f'{name}.json'], check=True, env=env, stdout=log, stderr=subprocess.STDOUT)
    with (output / 'consensus.log').open('w') as log:
        subprocess.run([sys.executable, root / 'scripts/consensus_human_family_pedigree.py',
            '--pair-predictions', output / 'pair_predictions.json', '--bonsai', output / 'fast.json',
            output / 'medium.json', '--topology-boost', '1', '--require-flat-degree-agreement',
            '--out-dir', output / 'consensus'], check=True, env=env, stdout=log, stderr=subprocess.STDOUT)
    with (output / 'consensus/candidate.csv').open() as handle:
        predictions = list(csv.DictReader(handle))
    verify_inputs(hashes)
    save_json(output / 'result.json', {'state': 'COMPLETED', 'pair': score_predictions(truth, predicted),
        'consensus': score_predictions(truth, predictions), 'submission_ready': False})



def neighbor_joining(distances: np.ndarray) -> dict:
    """Standard unrooted neighbor joining; no pedigree truth enters the fit."""
    n = len(distances)
    if distances.shape != (n, n) or n < 2 or not np.isfinite(distances).all():
        raise ValueError('Need a finite square distance matrix for at least two samples')
    if not np.allclose(distances, distances.T) or not np.allclose(np.diag(distances), 0):
        raise ValueError('Distances must be symmetric with zero diagonal')
    graph = defaultdict(dict)
    active = list(range(n))
    matrix = distances.copy()
    next_id = n
    while len(active) > 2:
        size = len(active)
        totals = matrix.sum(axis=1)
        q = (size - 2) * matrix - totals[:, None] - totals[None, :]
        np.fill_diagonal(q, np.inf)
        i, j = np.unravel_index(np.argmin(q), q.shape)
        if i > j:
            i, j = j, i
        limb_i = 0.5 * matrix[i, j] + (totals[i] - totals[j]) / (2 * (size - 2))
        limb_j = matrix[i, j] - limb_i
        node = next_id
        next_id += 1
        for index, length in [(i, limb_i), (j, limb_j)]:
            graph[node][active[index]] = graph[active[index]][node] = float(max(0.0, length))
        remaining = [k for k in range(size) if k not in (i, j)]
        joined = np.array([0.5 * (matrix[i, k] + matrix[j, k] - matrix[i, j]) for k in remaining])
        updated = np.zeros((size - 1, size - 1))
        updated[:-1, :-1] = matrix[np.ix_(remaining, remaining)]
        updated[-1, :-1] = updated[:-1, -1] = joined
        active = [active[k] for k in remaining] + [node]
        matrix = updated
    graph[active[0]][active[1]] = graph[active[1]][active[0]] = float(max(0.0, matrix[0, 1]))
    return dict(graph)


def integer_tree(graph: dict, samples: int) -> dict:
    """Round each inferred edge to meioses, retaining terminal edges of at least one."""
    representative = {node: node for node in graph}
    def find(node):
        while representative[node] != node:
            node = representative[node]
        return node
    edges = []
    for left, neighbors in graph.items():
        for right, length in neighbors.items():
            if left >= right:
                continue
            rounded = max(1 if left < samples or right < samples else 0, int(np.floor(length + 0.5)))
            if rounded == 0:
                representative[find(right)] = find(left)
            else:
                edges.append((left, right, rounded))
    output = defaultdict(set)
    next_id = max(graph) + 1
    for left, right, length in edges:
        left, right = find(left), find(right)
        previous = left
        for _ in range(length - 1):
            node = next_id
            next_id += 1
            output[previous].add(node)
            output[node].add(previous)
            previous = node
        output[previous].add(right)
        output[right].add(previous)
    if set(range(samples)) - set(output):
        raise ValueError('Integer tree lost a genotype leaf')
    return dict(output)


def rooted_labels(graph: dict, root: int, samples: list[str]) -> list[dict] | None:
    from ghist_relatedness.labels import accepted_label
    paths = {}
    stack = [(root, -1, [root])]
    while stack:
        node, parent, path = stack.pop()
        paths[node] = path
        stack.extend((child, node, path + [child]) for child in sorted(graph[node]) if child != parent)
    result = []
    for a, b in itertools.combinations(range(len(samples)), 2):
        pa, pb = paths[a], paths[b]
        common = 0
        for left, right in zip(pa, pb):
            if left != right:
                break
            common += 1
        up, down = len(pa) - common, len(pb) - common
        label = accepted_label((up, down, 2), 'human')
        if label is None:
            return None
        result.append({'indID1': samples[a], 'indID2': samples[b], 'relation': label})
    return result


def global_tree_predictions(features: list[dict], probabilities: list[dict], minimum_cm: int = 5) -> tuple[list[dict], dict]:
    """Rank only valid roots by a composite score, not a calibrated posterior."""
    samples = sorted({row[key] for row in features for key in ['indID1', 'indID2']})
    index = {sample: i for i, sample in enumerate(samples)}
    expected = set(itertools.combinations(samples, 2))
    pairs = {tuple(sorted((row['indID1'], row['indID2']))): row for row in features}
    probs = {tuple(sorted((row['indID1'], row['indID2']))): row['probabilities'] for row in probabilities}
    if len(pairs) != len(features) or set(pairs) != expected or set(probs) != expected:
        raise ValueError('Global inference requires every unique genotype roster pair')
    distances = np.zeros((len(samples), len(samples)))
    for (a, b), row in pairs.items():
        rel = float(row[f'relatedness_ge_{minimum_cm}'])
        if not np.isfinite(rel) or rel < 0:
            raise ValueError('Relatedness must be finite and nonnegative')
        distance = np.clip(1 - np.log2(max(rel, 2 ** -9)), 2, 10)
        distances[index[a], index[b]] = distances[index[b], index[a]] = distance
    raw_graph = neighbor_joining(distances)
    graph = integer_tree(raw_graph, len(samples))
    candidates = []
    for root in sorted(node for node in graph if node >= len(samples)):
        rows = rooted_labels(graph, root, samples)
        if rows is None:
            continue
        score = sum(np.log(max(float(probs[(row['indID1'], row['indID2'])][row['relation']]), 1e-12)) for row in rows)
        candidates.append((float(score), root, rows))
    candidates.sort(key=lambda item: (-item[0], item[1]))
    if not candidates:
        return [{key: row[key] for key in ['indID1', 'indID2', 'relation']} for row in probabilities], {
            'state': 'NO_VALID_ROOT_FALLBACK_TO_PAIR_MODEL', 'valid_roots': 0}
    best = candidates[0]
    alternatives = []
    for pair_index, row in enumerate(best[2]):
        labels = sorted({candidate[2][pair_index]['relation'] for candidate in candidates})
        if len(labels) > 1:
            alternatives.append({'indID1': row['indID1'], 'indID2': row['indID2'],
                                 'admissible_labels': labels})
    return best[2], {'state': 'COMPLETED', 'root': best[1], 'valid_roots': len(candidates),
        'root_score_margin': best[0] - candidates[1][0] if len(candidates) > 1 else None,
        'root_invariant_pairs': len(best[2]) - len(alternatives), 'root_ambiguous_pairs': len(alternatives),
        'root_ambiguity': alternatives, 'uncertainty_calibrated': False,
        'abstention_policy': 'Only root-invariant labels qualify as structurally resolved; ranked labels remain exploratory',
        'score_interpretation': 'composite pair log score; not an independent likelihood or calibrated posterior',
        'threshold_cm': minimum_cm, 'raw_edges': [[a, b, v] for a in raw_graph for b, v in raw_graph[a].items() if a < b],
        'integer_edges': [[a, b] for a in graph for b in sorted(graph[a]) if a < b]}


def run_tree(args) -> None:
    cohort, output = args.cohort.resolve(), args.out_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    with (cohort / 'diploid_state.csv').open() as handle:
        features = list(csv.DictReader(handle))
    probabilities = json.loads((cohort / 'diploid_state.json').read_text())
    truth = json.loads(args.truth.read_text())
    source_hashes = file_hashes(args.project_root.resolve())
    record = {'candidate': 'neighbor-joining-composite-root-v1', 'thresholds_cm': [5, 10],
        'source_hashes': source_hashes,
        'source_hash': sha256(Path(__file__)), 'input_hashes': {
            str(path): sha256(path) for path in [cohort / 'diploid_state.csv', cohort / 'diploid_state.json', args.truth]},
        'submission_ready': False, 'stage': 'DEVELOPMENT'}
    save_json(output / 'protocol.json', record)
    results = {}
    for threshold in [5, 10]:
        predictions, diagnostics = global_tree_predictions(features, probabilities, threshold)
        save_json(output / f'predictions_{threshold}.json', predictions)
        save_json(output / f'tree_{threshold}.json', diagnostics)
        results[str(threshold)] = {'score': score_predictions(truth, predictions), 'tree': diagnostics}
    verify_inputs(source_hashes)
    verify_inputs(record['input_hashes'])
    save_json(output / 'result.json', {'results': results, 'submission_ready': False})


def merge_spans(spans: list[tuple]) -> list[tuple]:
    merged = []
    for start, end in sorted(set(spans)):
        if end <= start:
            raise ValueError('Invalid interval')
        if merged and start <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def phase_pair_features(chromosomes: dict, genetic_map: dict, genome_cm: float) -> dict:
    """Chromosome-wise homolog scores with IBD2 masking and interval union."""
    raw, states = [], []
    side_minor = {threshold: [0.0, 0.0] for threshold in [5, 10]}
    side_total = {threshold: [0.0, 0.0] for threshold in [5, 10]}
    ibd2_mass = {threshold: 0.0 for threshold in [5, 10]}
    for chrom, original in chromosomes.items():
        positions, cm = genetic_map[chrom]
        by_edge = defaultdict(list)
        for start, end, hap1, hap2 in original:
            by_edge[(hap1, hap2)].append((start, end))
        intervals = []
        for (hap1, hap2), spans in by_edge.items():
            for start, end in merge_spans(spans):
                if start < positions[0] or end > positions[-1]:
                    raise ValueError('Phase feature interval exceeds genetic map')
                bounds = np.interp([start, end], positions, cm)
                length = float(bounds[1] - bounds[0])
                if length >= 2:
                    intervals.append((start, end, hap1, hap2))
                    raw.append((chrom, length))
        for start, end, full in coverage_intervals(intervals):
            bounds = np.interp([start, end], positions, cm)
            length = float(bounds[1] - bounds[0])
            states.extend([(chrom, length)] * (2 if full else 1))
        for threshold in [5, 10]:
            events = defaultdict(lambda: defaultdict(int))
            for start, end, h1, h2 in intervals:
                lo, hi = np.interp([start, end], positions, cm)
                if hi - lo < threshold:
                    continue
                events[float(lo)][(h1, h2)] += 1
                events[float(hi)][(h1, h2)] -= 1
            active, previous = defaultdict(int), None
            masses = [[0.0, 0.0], [0.0, 0.0]]
            for point in sorted(events):
                edges = {edge for edge, count in active.items() if count > 0}
                if previous is not None and edges:
                    length = point - previous
                    full = ({(1, 1), (2, 2)} <= edges or {(1, 2), (2, 1)} <= edges)
                    if full:
                        ibd2_mass[threshold] += length
                    else:
                        for side in [0, 1]:
                            for homolog in {edge[side] for edge in edges}:
                                masses[side][homolog - 1] += length
                for edge, count in events[point].items():
                    active[edge] += count
                previous = point
            for side in [0, 1]:
                side_minor[threshold][side] += min(masses[side])
                side_total[threshold][side] += sum(masses[side])
    features = {f'raw_{key}': value for key, value in summarize_raw_lengths(raw, genome_cm).items()}
    features.update({f'state_{key}': value for key, value in summarize_raw_lengths(states, genome_cm).items()})
    for threshold in [5, 10]:
        fractions = sorted(side_minor[threshold][side] / side_total[threshold][side]
            if side_total[threshold][side] else 0.0 for side in [0, 1])
        minor_cm = sorted(side_minor[threshold])
        features.update({f'phase_minor_low_{threshold}': fractions[0],
            f'phase_minor_high_{threshold}': fractions[1],
            f'phase_minor_cm_low_{threshold}': minor_cm[0], f'phase_minor_cm_high_{threshold}': minor_cm[1],
            f'ibd2_fraction_{threshold}': ibd2_mass[threshold] / genome_cm})
    return features


def phase_training_definition(copies: int) -> tuple[str, dict]:
    from ghist_relatedness.labels import HUMAN_LABELS
    lines, pairs = [], {}
    for (up, down, _), label in HUMAN_LABELS.items():
        lines.append(f'def {label} {copies} {down + 1}')
        if up == down:
            lines.append(f'{down + 1} 1')
        else:
            lines.append(f'{up + 1} 1 2 1n')
            for generation in range(up + 2, down + 1):
                lines.append(f'{generation} 0 1')
            lines.append(f'{down + 1} 1 1')
        for replicate in range(1, copies + 1):
            if up == down:
                samples = [f'{label}{replicate}_g{up + 1}-b1-i1', f'{label}{replicate}_g{down + 1}-b2-i1']
            else:
                samples = [f'{label}{replicate}_g{up + 1}-b2-i1', f'{label}{replicate}_g{down + 1}-b1-i1']
            pairs[tuple(sorted(samples))] = label
    return '\n'.join(lines) + '\n', pairs


def read_exact_phase(path: Path, expected: dict) -> dict:
    grouped = defaultdict(lambda: defaultdict(list))
    with path.open() as handle:
        for line in handle:
            fields = line.split()
            if len(fields) != 11:
                raise ValueError('Expected Ped-sim haplotype-resolved eleven-column output')
            a, h1, b, h2, chrom, start, end, state, _, _, _ = fields
            if state == 'HBD':
                continue
            if state != 'IBD1':
                raise ValueError('Unknown phased Ped-sim state')
            h1, h2 = int(h1) + 1, int(h2) + 1
            if a > b:
                a, b, h1, h2 = b, a, h2, h1
            if (a, b) not in expected:
                raise ValueError(f'Unexpected printed simulation pair {(a, b)}')
            if int(end) > int(start):
                grouped[(a, b)][chrom].append((int(start), int(end), h1, h2))
    return grouped


def train_phase_candidate(args) -> None:
    import joblib
    from sklearn.ensemble import ExtraTreesClassifier

    from ghist_relatedness.labels import HUMAN_LABELS
    output, assets = args.out_dir.resolve(), args.asset_root.resolve()
    output.mkdir(parents=True, exist_ok=False)
    protocol = json.loads(args.protocol.read_text())
    if protocol['candidate'] != 'phase-pair-forest-v1':
        raise ValueError('Unknown phased policy')
    pedsim = assets / 'external/ped-sim/ped-sim'
    simmap = assets / 'external/ped-sim-map/refined_mf.simmap'
    map_path = assets / 'external/ped-sim-map/refined_average.plink.map'
    interference = assets / 'external/ped-sim/interfere/nu_p_campbell.tsv'
    hashes = file_hashes(args.project_root.resolve())
    hashes.update({str(path): sha256(path) for path in [args.protocol.resolve(), pedsim, simmap, map_path, interference]})
    save_json(output / 'preregistered_design.json', {'protocol': protocol, 'hashes': hashes})
    genetic_map, genome_cm = read_map(map_path), genome_length(map_path)
    datasets = {}
    names = None
    for stage in ['training', 'exact_diagnostic']:
        spec = protocol[stage]
        definition, expected = phase_training_definition(spec['copies_per_label'])
        definition_path = output / f'{stage}.def'
        definition_path.write_text(definition)
        with (output / f'{stage}.log').open('w') as log:
            subprocess.run([str(pedsim), '-d', str(definition_path), '-m', str(simmap),
                '-o', str(output / stage), '--intf', str(interference), '--seed', str(spec['seed']),
                '--seg_haps'], check=True, stdout=log, stderr=subprocess.STDOUT)
        grouped = read_exact_phase(output / f'{stage}.seg', expected)
        rows = []
        for pair, label in expected.items():
            features = phase_pair_features(grouped.get(pair, {}), genetic_map, genome_cm)
            if names is None:
                names = sorted(features)
            rows.append({'indID1': pair[0], 'indID2': pair[1], 'relation': label, **features})
        save_json(output / f'{stage}_features.json', rows)
        datasets[stage] = rows
        print(stage, len(rows), flush=True)
    training, diagnostic = datasets['training'], datasets['exact_diagnostic']
    summary = {}
    for name, feature_names in [('phase', names), ('phase_free', [key for key in names if not key.startswith('phase_')])]:
        classifier = ExtraTreesClassifier(n_estimators=500, min_samples_leaf=5,
            class_weight='balanced', random_state=protocol['model_seed'], n_jobs=args.workers,
            max_features=1.0)
        classifier.fit([[row[key] for key in feature_names] for row in training], [row['relation'] for row in training])
        probabilities = classifier.predict_proba([[row[key] for key in feature_names] for row in diagnostic])
        predictions = [{'indID1': row['indID1'], 'indID2': row['indID2'],
            'relation': str(classifier.classes_[int(np.argmax(values))]),
            'probabilities': dict(zip(map(str, classifier.classes_), map(float, values), strict=True))}
            for row, values in zip(diagnostic, probabilities, strict=True)]
        bundle = {'model': classifier, 'feature_names': feature_names, 'protocol_sha256': sha256(args.protocol),
            'candidate': 'phase-pair-forest-v1', 'stage': 'DEVELOPMENT', 'variant': name, 'submission_ready': False}
        joblib.dump(bundle, output / f'{name}_model.joblib')
        save_json(output / f'{name}_exact_predictions.json', predictions)
        summary[name] = score_predictions(diagnostic, predictions)
    verify_inputs(hashes)
    save_json(output / 'result.json', {'state': 'COMPLETED', 'exact_diagnostic': summary,
        'labels': list(HUMAN_LABELS.values()), 'model_hashes': {name: sha256(output / f'{name}_model.joblib') for name in ['phase', 'phase_free']},
        'submission_ready': False, 'submission_count': 0})


def predict_phase_candidate(args) -> None:
    import joblib
    output, assets, cohort = args.out_dir.resolve(), args.asset_root.resolve(), args.cohort.resolve()
    output.mkdir(parents=True, exist_ok=False)
    model = joblib.load(args.model)
    map_path = assets / 'external/ped-sim-map/refined_average.plink.map'
    genetic_map, genome_cm = read_map(map_path), genome_length(map_path)
    roster = genotype_roster(cohort / 'family.normalized.vcf.gz')
    grouped, _ = read_haplotype_segments(cohort / 'detected.ibd.gz')
    by_pair = defaultdict(dict)
    for (a, b, chrom), intervals in grouped.items():
        by_pair[(a, b)][chrom] = intervals
    hashes = file_hashes(args.project_root.resolve())
    hashes.update({str(path): sha256(path) for path in [args.model.resolve(), map_path,
        cohort / 'detected.ibd.gz', cohort / 'family.normalized.vcf.gz', cohort / 'truth.json']})
    save_json(output / 'protocol.json', {'candidate': model['candidate'], 'variant': model['variant'],
        'hashes': hashes, 'stage': 'DEVELOPMENT'})
    features = []
    for a, b in itertools.combinations(sorted(roster), 2):
        values = phase_pair_features(by_pair.get((a, b), {}), genetic_map, genome_cm)
        features.append({'indID1': a, 'indID2': b, **values})
    probabilities = model['model'].predict_proba([[row[key] for key in model['feature_names']] for row in features])
    predictions = [{'indID1': row['indID1'], 'indID2': row['indID2'],
        'relation': str(model['model'].classes_[int(np.argmax(values))]),
        'probabilities': dict(zip(map(str, model['model'].classes_), map(float, values), strict=True))}
        for row, values in zip(features, probabilities, strict=True)]
    save_json(output / 'features.json', features)
    save_json(output / 'predictions.json', predictions)
    truth = json.loads((cohort / 'truth.json').read_text())
    verify_inputs(hashes)
    save_json(output / 'result.json', {'state': 'COMPLETED', 'score': score_predictions(truth, predictions),
        'submission_ready': False, 'submission_count': 0})

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['diagnose', 'whole-family', 'cohort', 'tree', 'phase-train', 'phase-predict'])
    parser.add_argument('--cohort', type=Path)
    parser.add_argument('--truth', type=Path)
    parser.add_argument('--model', type=Path)
    parser.add_argument('--asset-root', type=Path, required=True)
    parser.add_argument('--out-dir', type=Path, required=True)
    parser.add_argument('--project-root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--protocol', type=Path)
    parser.add_argument('--index', type=int)
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):
        parser.error('Execute model and full-genome work in scheduled Slurm jobs')
    if args.command == 'diagnose':
        print(json.dumps(evaluate_representations(args.cohort, args.asset_root, args.out_dir), indent=2))
    elif args.command == 'phase-train':
        train_phase_candidate(args)
    elif args.command == 'phase-predict':
        predict_phase_candidate(args)
    elif args.command == 'tree':
        if args.truth is None:
            parser.error('Tree development scoring requires truth separate from inference')
        run_tree(args)
    elif args.command == 'whole-family':
        run_whole_family(args)
    else:
        if args.protocol is None or args.index is None:
            parser.error('Confirmation requires a preregistered protocol and cohort index')
        run_cohort(args)


if __name__ == '__main__':
    main()
