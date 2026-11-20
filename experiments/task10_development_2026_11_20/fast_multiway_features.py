"""Exact feature optimization candidate; independent of the frozen v2 implementation.

Three-group sharing requires simultaneous edges of at least two group-pair
types. Outside those intervals, every connected component lies in at most two
groups, so its conflicts are covered by the union of precomputed pair conflicts.
No biological features, thresholds, groups or model parameters are changed.
"""
from __future__ import annotations

import bisect
import itertools
from collections import defaultdict

import numpy as np
from multiway_candidate import inferred_groups

from scripts.develop_task10_detected_ibd import merge_spans, phase_pair_features


NAMES = ['union', 'copies', 'multiple', 'conflict', 'complete']


def union_cm(spans):
    merged = []
    for start, end in sorted(set(spans)):
        if end <= start:
            raise ValueError('Invalid genetic interval')
        if merged and start <= merged[-1][1]:
            merged[-1] = merged[-1][0], max(merged[-1][1], end)
        else:
            merged.append((start, end))
    return merged


def overlap_two_types(interval_sets):
    events = defaultdict(int)
    for intervals in interval_sets:
        for start, end in union_cm(intervals):
            events[start] += 1
            events[end] -= 1
    count, previous, result = 0, None, []
    for position in sorted(events):
        if previous is not None and count >= 2 and position > previous:
            result.append((previous, position))
        count += events[position]
        previous = position
    return union_cm(result)


def sweep(edges, membership, group_count, windows=None):
    events = defaultdict(lambda: defaultdict(int))
    stops = [end for _, end in windows] if windows is not None else None
    for left, right, start, end in edges:
        intervals = [(start, end)]
        if windows is not None:
            intervals = []
            index = bisect.bisect_right(stops, start)
            while index < len(windows) and windows[index][0] < end:
                lo, hi = max(start, windows[index][0]), min(end, windows[index][1])
                if hi > lo:
                    intervals.append((lo, hi))
                index += 1
        for lo, hi in intervals:
            events[lo][(left, right)] += 1
            events[hi][(left, right)] -= 1
    totals, conflict_spans = dict.fromkeys(NAMES, 0.), []
    active, previous = defaultdict(int), None
    required = set(itertools.combinations(sorted(set(membership.values())), 2))
    for position in sorted(events):
        current = [edge for edge, count in active.items() if count > 0]
        if previous is not None and current:
            parents = {node: node for edge in current for node in edge}

            def find(node):
                while parents[node] != node:
                    parents[node] = parents[parents[node]]
                    node = parents[node]
                return node

            for left, right in current:
                parents[find(left)] = find(right)
            components, types = defaultdict(set), defaultdict(set)
            for node in parents:
                components[find(node)].add(node)
            for left, right in current:
                a, b = membership[left // 2], membership[right // 2]
                if a != b:
                    types[find(left)].add(tuple(sorted([a, b])))
            shared, complete, conflict = 0, 0, False
            for root, nodes in components.items():
                people = [node // 2 for node in nodes]
                if len(people) != len(set(people)):
                    conflict = True
                    continue
                if len({membership[person] for person in people}) == group_count:
                    shared += 1
                    complete += required <= types[root]
            length = position - previous
            totals['union'] += length * (shared > 0)
            totals['copies'] += length * shared
            totals['multiple'] += length * (shared >= 2)
            totals['conflict'] += length * conflict
            totals['complete'] += length * (complete > 0)
            if conflict:
                conflict_spans.append((previous, position))
        for edge, count in events[position].items():
            active[edge] += count
        previous = position
    return totals, union_cm(conflict_spans)


class ComponentCache:
    def __init__(self, groups, segments, genetic_map, genome_cm):
        self.groups, self.genome_cm = groups, genome_cm
        samples = sorted(sample for group in groups for sample in group)
        identifiers = {sample: index for index, sample in enumerate(samples)}
        self.membership = {identifiers[sample]: index for index, group in enumerate(groups) for sample in group}
        self.edges = {threshold: defaultdict(lambda: defaultdict(list)) for threshold in [5, 10]}
        self.cross = {threshold: defaultdict(dict) for threshold in [5, 10]}
        self.pairs, self.triples = {}, {}
        for (a, b), chromosomes in segments.items():
            if a not in identifiers or b not in identifiers:
                continue
            group_pair = tuple(sorted([self.membership[identifiers[a]], self.membership[identifiers[b]]]))
            for chrom, rows in chromosomes.items():
                by_edge = defaultdict(list)
                for start, end, h1, h2 in rows:
                    by_edge[(2 * identifiers[a] + h1 - 1, 2 * identifiers[b] + h2 - 1)].append((start, end))
                for edge, spans in by_edge.items():
                    for start, end in merge_spans(spans):
                        lo, hi = map(float, np.interp([start, end], *genetic_map[chrom]))
                        for threshold in [5, 10]:
                            if hi - lo >= threshold:
                                self.edges[threshold][group_pair][chrom].append((*edge, lo, hi))
        for threshold, groups_edges in self.edges.items():
            for pair, chromosomes in groups_edges.items():
                if pair[0] != pair[1]:
                    self.cross[threshold][pair] = {chrom: union_cm([(start, end) for _, _, start, end in edges])
                        for chrom, edges in chromosomes.items()}

    def relevant_edges(self, indices, threshold):
        result = defaultdict(list)
        for pair in itertools.combinations_with_replacement(indices, 2):
            for chrom, edges in self.edges[threshold].get(pair, {}).items():
                result[chrom].extend(edges)
        return result

    def integral(self, indices, threshold):
        indices = tuple(sorted(indices))
        key = indices, threshold
        cache = self.pairs if len(indices) == 2 else self.triples
        if key in cache:
            return cache[key][0]
        if len(indices) not in [2, 3] or len(set(indices)) != len(indices):
            raise ValueError('Expected two or three distinct inferred groups')
        membership = {person: group for person, group in self.membership.items() if group in indices}
        totals, conflicts = dict.fromkeys(NAMES, 0.), {}
        if len(indices) == 2:
            edges = self.relevant_edges(indices, threshold)
            for chrom, rows in edges.items():
                values, spans = sweep(rows, membership, 2)
                conflicts[chrom] = spans
                for name, value in values.items():
                    totals[name] += value
        else:
            pairs = list(itertools.combinations(indices, 2))
            for pair in pairs:
                self.integral(pair, threshold)
            edge_sources = [self.edges[threshold].get(pair, {})
                for pair in itertools.combinations_with_replacement(indices, 2)]
            chromosomes = sorted({chrom for source in edge_sources for chrom in source})
            for chrom in chromosomes:
                windows = overlap_two_types([self.cross[threshold].get(pair, {}).get(chrom, []) for pair in pairs])
                rows = [edge for source in edge_sources for edge in source.get(chrom, [])] if windows else []
                values, new_conflicts = sweep(rows, membership, 3, windows) if windows else (dict.fromkeys(NAMES, 0.), [])
                existing = [span for pair in pairs for span in self.pairs[(pair, threshold)][1].get(chrom, [])]
                conflicts[chrom] = union_cm(existing + new_conflicts)
                for name, value in values.items():
                    if name != 'conflict':
                        totals[name] += value
                totals['conflict'] += sum(end - start for start, end in conflicts[chrom])
        result = {name: value / self.genome_cm for name, value in totals.items()}
        cache[key] = result, conflicts if len(indices) == 2 else {}
        return result


def family_features(samples, segments, genetic_map, genome_cm, phase_model):
    rows = [{'indID1': a, 'indID2': b, **phase_pair_features(segments.get((a, b), {}), genetic_map, genome_cm)}
        for a, b in itertools.combinations(sorted(samples), 2)]
    groups, frozen_probabilities = inferred_groups(samples, rows, phase_model)
    membership = {sample: index for index, group in enumerate(groups) for sample in group}
    engine = ComponentCache(groups, segments, genetic_map, genome_cm)
    for row in rows:
        a, b = membership[row['indID1']], membership[row['indID2']]
        low, high = sorted([len(groups[a]), len(groups[b])])
        row.update(context_same_group=float(a == b), context_group_size_low=float(low),
            context_group_size_high=float(high), context_groups=float(len(groups)))
        for threshold in [5, 10]:
            pooled = engine.integral([a, b], threshold) if a != b else dict.fromkeys(NAMES, 0.)
            for name, value in pooled.items():
                row[f'context_pool_{name}_{threshold}'] = value
            triples = [engine.integral([a, b, c], threshold) for c in range(len(groups)) if c not in [a, b]] if a != b else []
            for name in [*NAMES, 'relative_union']:
                values = [triple[name] if name != 'relative_union' else triple['union'] / max(pooled['union'], 1e-12)
                    for triple in triples]
                statistics = [float(np.mean(values)), *map(float, np.quantile(values, [0., .25, .5, .75, 1.]))] if values else [0.] * 6
                for statistic, value in zip(['mean', 'min', 'q25', 'median', 'q75', 'max'], statistics, strict=True):
                    row[f'context_triple_{name}_{statistic}_{threshold}'] = value
    return rows, groups, frozen_probabilities
