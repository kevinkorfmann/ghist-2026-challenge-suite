"""Independent pedigree checks for the frozen sparse-context stress design."""
import json
from collections import deque
from itertools import combinations, product
from pathlib import Path

import numpy as np
import pytest

from scripts import develop_task11_pedigree as dev


def ancestor_depths(parents, sample):
    queue = deque([(sample, 0)])
    depths = {}
    while queue:
        node, depth = queue.popleft()
        if node in depths and depths[node] <= depth:
            continue
        depths[node] = depth
        if parents[node] is not None:
            queue.extend((parent, depth + 1) for parent in parents[node])
    return depths


@pytest.mark.parametrize("label", dev.LABELS)
@pytest.mark.parametrize("sizes", list(product(range(1, 5), repeat=2)))
@pytest.mark.parametrize("third,complex_topology", list(product([False, True], repeat=2)))
def test_every_sparse_pair_has_independent_ancestor_derived_label(label, sizes, third, complex_topology):
    pedigree = dev.focal_pedigree(label, sizes, third, complex_topology)
    pairs = pedigree["pairs"]
    assert {frozenset((a, b)) for a, b, _ in pairs} == {
        frozenset(pair) for pair in combinations(pedigree["observed"], 2)
    }
    assert len(pairs) == len({frozenset((a, b)) for a, b, _ in pairs})
    assert label in {relation for _, _, relation in pairs}
    for a, b, relation in pairs:
        left = ancestor_depths(pedigree["parents"], a)
        right = ancestor_depths(pedigree["parents"], b)
        assert a not in right and b not in left
        common = left.keys() & right.keys()
        nearest_sum = min(left[node] + right[node] for node in common)
        nearest = [node for node in common if left[node] + right[node] == nearest_sum]
        observed_depths = {tuple(sorted((left[node], right[node]))) for node in nearest}
        d1, d2, shared = dev.DEPTHS[relation]
        assert observed_depths == {(d1, d2)}
        assert len(nearest) == shared


def test_stress_seeds_are_disjoint_and_actual_topologies_include_singletons():
    design = json.loads((Path(__file__).parent / "sparse_context_stress_design_v1.json").read_text())
    seeds = [cohort["seed"] for cohort in design["cohorts"]]
    assert len(seeds) == len(set(seeds)) == 8
    assert not set(seeds) & set(design["existing_and_reserved_seed_union"])
    assert design["gates"] == dev.GATES
    for cohort in design["cohorts"]:
        rng = np.random.default_rng(cohort["seed"])
        sizes = rng.integers(1, 5, (8, 2))
        assert 1 in sizes
        assert 4 in sizes
    assert sum(cohort["complex_topology"] for cohort in design["cohorts"]) == 4
