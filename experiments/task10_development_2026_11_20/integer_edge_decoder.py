"""Fixed bounded integer-edge optimization on an inferred NJ topology only."""
from __future__ import annotations

import argparse
import itertools
import json
import math
import os
from collections import defaultdict, deque
from pathlib import Path

import numpy as np
from calibrate_phase import metrics
from ghist_relatedness.labels import HUMAN_LABELS

from scripts.develop_task10_detected_ibd import rooted_labels, save_json, sha256, verify_inputs


def core_tree(integer_edges, sample_count):
    graph = defaultdict(set)
    for a, b in integer_edges:
        if a == b or b in graph[a]:
            raise ValueError('Invalid tree edge')
        graph[a].add(b)
        graph[b].add(a)
    leaves = {node for node, neighbors in graph.items() if len(neighbors) == 1}
    if len(integer_edges) != len(graph) - 1 or leaves != set(range(sample_count)):
        raise ValueError('Expected a tree with every sample as a leaf')
    seen, queue = set(), [0]
    while queue:
        node = queue.pop()
        if node not in seen:
            seen.add(node)
            queue.extend(graph[node] - seen)
    if seen != set(graph):
        raise ValueError('Disconnected tree')
    core = {node for node in graph if node < sample_count or len(graph[node]) != 2}
    edges = []
    for a in sorted(core):
        for next_node in sorted(graph[a]):
            previous, current, length = a, next_node, 1
            while current not in core:
                following = next(node for node in graph[current] if node != previous)
                previous, current = current, following
                length += 1
            if a < current:
                edges.append((a, current, length))
    return sorted(edges)


def expand_tree(edges, lengths):
    graph = defaultdict(set)
    next_node = 1 + max(node for a, b, _ in edges for node in [a, b])
    for (a, b, _), length in zip(edges, lengths, strict=True):
        previous = a
        for _ in range(int(length) - 1):
            graph[previous].add(next_node)
            graph[next_node].add(previous)
            previous = next_node
            next_node += 1
        graph[previous].add(b)
        graph[b].add(previous)
    return dict(graph)


def pair_paths(edges, sample_count):
    graph = defaultdict(list)
    for index, (a, b, _) in enumerate(edges):
        graph[a].append((b, index))
        graph[b].append((a, index))
    pairs = list(itertools.combinations(range(sample_count), 2))
    matrix = np.zeros((len(pairs), len(edges)), dtype=np.int16)
    for row, (a, b) in enumerate(pairs):
        queue, paths = deque([a]), {a: []}
        while b not in paths:
            node = queue.popleft()
            for neighbor, edge in graph[node]:
                if neighbor not in paths:
                    paths[neighbor] = paths[node] + [edge]
                    queue.append(neighbor)
        matrix[row, paths[b]] = 1
    return matrix


def valid_root_rows(graph, samples):
    """Prune only roots violating the accepted-label maximum of five meioses.

    Every internal root has at least two descendant branches containing leaves;
    hence a leaf farther than five gives a pair whose LCA is this root and whose
    label is unsupported. This is a label-domain constraint, not a root prior.
    """
    for root in sorted(node for node in graph if node >= len(samples)):
        distances, queue = {root: 0}, deque([root])
        while queue:
            node = queue.popleft()
            for neighbor in graph[node]:
                if neighbor not in distances:
                    distances[neighbor] = distances[node] + 1
                    queue.append(neighbor)
        if max(distances[index] for index in range(len(samples))) > 5:
            continue
        rows = rooted_labels(graph, root, samples)
        if rows is not None:
            yield root, rows


def optimize(predictions, integer_edges, max_moves=20, max_feasibility_checks=200):
    samples = sorted({row[key] for row in predictions for key in ['indID1', 'indID2']})
    keyed = {tuple(sorted([row['indID1'], row['indID2']])): row for row in predictions}
    pairs = list(itertools.combinations(samples, 2))
    if set(keyed) != set(pairs) or len(keyed) != len(predictions):
        raise ValueError('Expected every unique roster pair')
    meioses = {label: up + down for (up, down, _), label in HUMAN_LABELS.items()}
    probabilities = np.zeros((len(pairs), 11))
    for index, pair in enumerate(pairs):
        values = keyed[pair]['probabilities']
        if set(values) != set(meioses) or any(not math.isfinite(p) or p < 0 for p in values.values()):
            raise ValueError('Invalid label probabilities')
        if not math.isclose(sum(values.values()), 1., abs_tol=1e-10, rel_tol=1e-10):
            raise ValueError('Probabilities must sum to one')
        for label, value in values.items():
            probabilities[index, meioses[label]] += value
    logs = np.log(np.maximum(probabilities, 1e-12))
    edges = core_tree(integer_edges, len(samples))
    paths = pair_paths(edges, len(samples))
    lengths = np.asarray([length for _, _, length in edges], dtype=np.int16)
    distance = paths @ lengths
    if np.any(distance < 2) or np.any(distance > 10):
        raise ValueError('Initial inferred tree lies outside accepted meiosis counts')
    current = math.fsum(logs[index, degree] for index, degree in enumerate(distance))
    initial, moves, feasibility = current, [], 0
    visited = {tuple(lengths)}
    stop = 'MOVE_CAP'
    for _ in range(max_moves):
        candidates = []
        for edge, length in enumerate(lengths):
            affected = np.flatnonzero(paths[:, edge])
            old = distance[affected]
            for change in [-1, 1]:
                if not 1 <= length + change <= 10:
                    continue
                changed = old + change
                if np.any(changed < 2) or np.any(changed > 10):
                    continue
                score = current + math.fsum(logs[index, new] - logs[index, prior]
                    for index, prior, new in zip(affected, old, changed, strict=True))
                if score > current + 1e-12:
                    candidates.append((score, edge, change))
        candidates.sort(key=lambda item: (-item[0], item[1], item[2]))
        accepted = False
        for score, edge, change in candidates:
            if feasibility >= max_feasibility_checks:
                stop = 'FEASIBILITY_CAP'
                break
            proposal = lengths.copy()
            proposal[edge] += change
            if tuple(proposal) in visited:
                continue
            feasibility += 1
            graph = expand_tree(edges, proposal)
            if next(valid_root_rows(graph, samples), None) is None:
                continue
            old_length = int(lengths[edge])
            lengths, distance, current = proposal, paths @ proposal, score
            visited.add(tuple(lengths))
            moves.append({'edge_endpoints': list(edges[edge][:2]), 'old_length': old_length,
                'new_length': int(proposal[edge]), 'composite_degree_log_score': current})
            accepted = True
            break
        if not accepted:
            if stop != 'FEASIBILITY_CAP':
                stop = 'NO_FEASIBLE_IMPROVING_SINGLE_EDGE_MOVE'
            break
    graph = expand_tree(edges, lengths)
    roots = []
    for root, rows in valid_root_rows(graph, samples):
        score = sum(np.log(max(float(keyed[(row['indID1'], row['indID2'])]['probabilities'][row['relation']]), 1e-12))
            for row in rows)
        roots.append((float(score), root, rows))
    roots.sort(key=lambda item: (-item[0], item[1]))
    if not roots:
        raise ValueError('Optimization lost all valid roots')
    best = roots[0]
    ambiguous = [{'indID1': row['indID1'], 'indID2': row['indID2'],
        'admissible_labels': sorted({candidate[2][index]['relation'] for candidate in roots})}
        for index, row in enumerate(best[2])
        if len({candidate[2][index]['relation'] for candidate in roots}) > 1]
    return best[2], {'state': 'COMPLETED', 'root': best[1], 'valid_roots': len(roots),
        'initial_degree_log_score': initial, 'optimized_degree_log_score': current,
        'moves': moves, 'feasibility_checks': feasibility, 'stop': stop,
        'root_score_margin': best[0] - roots[1][0] if len(roots) > 1 else None,
        'root_invariant_pairs': len(best[2]) - len(ambiguous), 'root_ambiguous_pairs': len(ambiguous),
        'root_ambiguity': ambiguous,
        'integer_edges': [[a, b] for a in graph for b in sorted(graph[a]) if a < b],
        'topology_fixed': True, 'uncertainty_calibrated': False,
        'score_interpretation': 'Composite degree log score for edge lengths, full-label composite score for roots; neither is an independent likelihood or calibrated posterior',
        'uncertainty_scope': 'Root ambiguity is conditional on selected edge lengths and does not cover other edge-length assignments or other topologies'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', type=Path, required=True)
    parser.add_argument('--protocol', type=Path, required=True)
    parser.add_argument('--panel-root', type=Path, required=True)
    parser.add_argument('--initial-tree-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):
        parser.error('Use Slurm')
    protocol = json.loads(args.protocol.read_text())
    if sha256(Path(__file__)) != protocol['decoder_source_sha256']:
        raise ValueError('Decoder source changed')
    if sha256(args.source_root / 'source_manifest.json') != protocol['candidate_source_manifest_sha256']:
        raise ValueError('Frozen source changed')
    verify_inputs(json.loads((args.source_root / 'source_manifest.json').read_text()))
    evaluation = args.panel_root / 'model_evaluation'
    previous = json.loads((evaluation / 'result.json').read_text())
    bound = {str(path): sha256(path) for path in [args.protocol, Path(__file__),
        evaluation / 'result.json', evaluation / 'model_calibration_freeze.json', args.initial_tree_root / 'result.json']}
    args.output.mkdir(parents=True, exist_ok=False)
    scores, gates = {}, {}
    for variant in protocol['variants']:
        source_path = evaluation / f'{variant}_development_scorer_input.json'
        tree_path = args.initial_tree_root / f'{variant}_trees.json'
        for path in [source_path, tree_path]:
            bound[str(path)] = sha256(path)
        source = json.loads(source_path.read_text())
        initial = {row['cohort_id']: row for row in json.loads(tree_path.read_text())}
        cohorts, trees = [], []
        for cohort in source['cohorts']:
            predictions, tree = optimize(cohort['predictions'], initial[cohort['cohort_id']]['integer_edges'],
                protocol['max_accepted_moves'], protocol['max_feasibility_checks'])
            cohorts.append({'cohort_id': cohort['cohort_id'], 'truth': cohort['truth'], 'predictions': predictions})
            trees.append({'cohort_id': cohort['cohort_id'], **tree})
        save_json(args.output / f'{variant}_scorer_input.json', {'labels': source['labels'], 'cohorts': cohorts})
        save_json(args.output / f'{variant}_predictions.json', [{**row, 'cohort_id': cohort['cohort_id']}
            for cohort in cohorts for row in cohort['predictions']])
        save_json(args.output / f'{variant}_trees.json', trees)
        score = scores[variant] = metrics(cohorts, source['labels'])
        baseline = previous['baseline']
        gates[variant] = {'all_labels': score['complete_label_coverage'], 'exact_accuracy': score['exact_accuracy'] >= .8,
            'macro_recall': score['macro_recall'] >= .7, 'minimum_label_recall': score['minimum_label_recall'] >= .5,
            'gain_over_baseline': score['exact_accuracy'] - baseline['exact_accuracy'] >= .1,
            'no_cohort_regression': all(value >= baseline['cohort_accuracy'][cohort]
                for cohort, value in score['cohort_accuracy'].items())}
    verify_inputs(bound)
    save_json(args.output / 'result.json', {'stage': 'REUSED_DEVELOPMENT_DIAGNOSTIC',
        'scores': scores, 'gates': gates, 'input_sha256': bound, 'protocol_sha256': sha256(args.protocol),
        'original_v2_policy_changed': False, 'uncertainty_calibrated': False,
        'submission_ready': False, 'advance_to_confirmation': False})


if __name__ == '__main__':
    main()
