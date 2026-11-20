"""Post-outcome oracle diagnosis only; never predicts or selects a usable root."""
import argparse
import itertools
import json
from collections import Counter, deque
from pathlib import Path

from ghist_relatedness.labels import HUMAN_LABELS
from scripts.develop_task10_detected_ibd import rooted_labels, save_json, sha256


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    source = args.evidence / 'context_scorer_input.json'
    tree_path = args.evidence / 'context_trees.json'
    cohorts = {row['cohort_id']: row for row in json.loads(source.read_text())['cohorts']}
    relation = {label: up + down for (up, down, _), label in HUMAN_LABELS.items()}
    reports = []
    for tree in json.loads(tree_path.read_text()):
        cohort = cohorts[tree['cohort_id']]
        truth = {(row['indID1'], row['indID2']): row['relation'] for row in cohort['truth']}
        samples = sorted({sample for pair in truth for sample in pair})
        graph = {}
        for a, b in tree['integer_edges']:
            graph.setdefault(a, set()).add(b)
            graph.setdefault(b, set()).add(a)
        roots = []
        for root in sorted(node for node in graph if node >= len(samples)):
            rows = rooted_labels(graph, root, samples)
            if rows is not None:
                roots.append(sum(row['relation'] == truth[(row['indID1'], row['indID2'])] for row in rows))
        distances = {}
        for a in range(len(samples)):
            distances[a], queue = {a: 0}, deque([a])
            while queue:
                node = queue.popleft()
                for neighbor in graph[node]:
                    if neighbor not in distances[a]:
                        distances[a][neighbor] = distances[a][node] + 1
                        queue.append(neighbor)
        errors = Counter()
        for a, b in itertools.combinations(range(len(samples)), 2):
            label = truth[(samples[a], samples[b])]
            errors[(relation[label], distances[a][b])] += 1
        reports.append({'cohort_id': tree['cohort_id'], 'pairs': len(truth),
            'predicted_correct': sum(row['relation'] == truth[(row['indID1'], row['indID2'])] for row in cohort['predictions']),
            'oracle_max_correct_among_inferred_valid_roots': max(roots),
            'inferred_meiosis_distance_correct': sum(count for (actual, inferred), count in errors.items() if actual == inferred),
            'distance_confusion': [{'true_meioses': a, 'inferred_meioses': b, 'count': count} for (a, b), count in sorted(errors.items())]})
    save_json(args.output, {'stage': 'POST_OUTCOME_ORACLE_DIAGNOSTIC_NOT_PREDICTIONS',
        'cohorts': reports, 'totals': {key: sum(row[key] for row in reports) for key in ['pairs', 'predicted_correct',
            'oracle_max_correct_among_inferred_valid_roots', 'inferred_meiosis_distance_correct']},
        'input_sha256': {str(path): sha256(path) for path in [source, tree_path]},
        'source_sha256': sha256(Path(__file__)), 'submission_ready': False})
    print(json.dumps(reports, indent=2))


if __name__ == '__main__':
    main()
