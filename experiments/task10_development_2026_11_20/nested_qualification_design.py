"""Prospective, outcome-independent nested pedigree qualification definitions."""
from __future__ import annotations

import itertools
from collections import Counter, defaultdict

from ghist_relatedness.labels import HUMAN_LABELS


def nested_shapes():
    base = [(clade, *tail) for clade in [0, 1]
        for tail in [(0,), (1, 0), (1, 1, 0), (1, 1, 1, 0)]]
    return [base,
        base + [(0, 2, 0), (1, 2, 0)],
        base + [(0, 1, 2, 0), (1, 1, 2, 0)],
        base + [(0, 2, 0), (1, 2, 0), (0, 1, 2, 0), (1, 1, 2, 0)]]


def large_nested_shape():
    leaves = list(nested_shapes()[0])
    for branch in range(2, 109):
        depth = 2 + (branch - 2) % 4
        leaves.append((branch, *([0] * (depth - 1))))
    if len(leaves) != 115:
        raise ValueError('Expected 115 terminal branches')
    return leaves


def compile_definition(leaves, name, people=None):
    leaves = sorted(map(tuple, leaves))
    if len(set(leaves)) != len(leaves) or any(not 2 <= len(path) <= 5 for path in leaves):
        raise ValueError('Expected unique terminal paths of depths two through five')
    if any(a[:len(b)] == b for a, b in itertools.permutations(leaves, 2)):
        raise ValueError('Observed terminal paths cannot be ancestors of each other')
    nodes = {depth: sorted({path[:depth] for path in leaves if len(path) >= depth})
        for depth in range(1, max(map(len, leaves)) + 1)}
    index = {depth: {path: branch for branch, path in enumerate(paths, 1)} for depth, paths in nodes.items()}
    lines = [f'def {name} 1 {max(map(len, leaves)) + 1}']
    for depth, paths in nodes.items():
        parent_groups = defaultdict(list)
        for path in paths:
            parent_groups[1 if depth == 1 else index[depth - 1][path[:-1]]].append(index[depth][path])
        assignments = [','.join(map(str, branches)) + f':{parent}' for parent, branches in parent_groups.items()]
        omitted = [f'{index[depth][path]}n' for path in paths if path not in leaves]
        lines.append(f'{depth + 1} 2 {len(paths)} ' + ' '.join(assignments + omitted))
    ordered = [(f'{name}1_g{len(path) + 1}-b{index[len(path)][path]}-i{child}', path)
        for path in leaves for child in [1, 2]]
    if people is not None:
        if people < 2 or people > len(ordered):
            raise ValueError('Invalid declared genotype roster size')
        ordered = ordered[:people]
    expected = {}
    for (a, left), (b, right) in itertools.combinations(ordered, 2):
        if left == right:
            relation = (1, 1, 2)
        else:
            shared = next((i for i, (x, y) in enumerate(zip(left, right, strict=False)) if x != y), min(len(left), len(right)))
            relation = (*sorted([len(left) - shared, len(right) - shared]), 2)
        expected[tuple(sorted([a, b]))] = HUMAN_LABELS[relation]
    return '\n'.join(lines) + '\n', [person for person, _ in ordered], expected


def design_summary(leaves, name, people=None):
    definition, roster, expected = compile_definition(leaves, name, people)
    return {'name': name, 'terminal_paths': leaves, 'people': len(roster), 'pairs': len(expected),
        'label_pair_counts': dict(sorted(Counter(expected.values()).items())),
        'definition': definition, 'declared_roster': roster}
