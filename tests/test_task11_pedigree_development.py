from itertools import combinations

import numpy as np
import pytest

from scripts import develop_task11_pedigree as dev


def distances(parents, sample):
    found = {}
    queue = [(sample, 0)]
    while queue:
        node, distance = queue.pop(0)
        if node in found and found[node] <= distance:
            continue
        found[node] = distance
        if parents[node] is not None:
            queue.extend((p, distance + 1) for p in parents[node])
    return found


@pytest.mark.parametrize("label", dev.LABELS)
@pytest.mark.parametrize("third", [False, True])
def test_focal_truth_matches_parent_dag_and_covers_every_observed_pair(label, third):
    p = dev.focal_pedigree(label, (2, 3), third)
    observed = p["observed"]
    assert len(p["pairs"]) == len(list(combinations(observed, 2)))
    assert len({tuple(sorted((a, b))) for a, b, _ in p["pairs"]}) == len(p["pairs"])
    for a, b, relation in p["pairs"]:
        da, db = distances(p["parents"], a), distances(p["parents"], b)
        assert a not in db and b not in da
        common = set(da) & set(db)
        shortest = min(da[n] + db[n] for n in common)
        nearest = [n for n in common if da[n] + db[n] == shortest]
        depths = tuple(sorted((da[nearest[0]], db[nearest[0]])))
        left, right, n = dev.DEPTHS[relation]
        assert depths == (left, right)
        assert len(nearest) == n


@pytest.mark.parametrize("kind", ["unrelated", "first_cousins", "siblings"])
def test_background_inbreeding_keeps_topological_order_and_observed_ids(kind):
    p = dev.focal_pedigree("avuncular_once_removed", (3, 4), True)
    parents = dev.expand_inbreeding(p, kind)
    seen = set()
    for name, pair in parents.items():
        if pair is not None:
            assert set(pair) <= seen
        seen.add(name)
    assert set(p["observed"]) <= seen
    a, b = p["central_founders"]
    common = set(distances(parents, a)) & set(distances(parents, b))
    assert bool(common) == (kind != "unrelated")


def row(a, b, relation="siblings", support=1):
    return {"indID1": a, "indID2": b, "relation": relation, "support": support, "family": 0}


def test_unconditional_score_retains_missing_and_extra_pairs():
    truth = {"pairs": [row("a", "b"), row("a", "c", "half_sibling")]}
    result = dev.score(truth, [row("a", "b"), row("x", "y")])
    assert result["accuracy"] == 0.5
    assert result["missing_true_pairs"] == 1
    assert result["extra_cross_family_pairs"] == 1
    assert result["recall"]["half_sibling"] == 0
    with pytest.raises(ValueError, match="Duplicate"):
        dev.score(truth, [row("a", "b"), row("b", "a")])


def test_temperature_keeps_ranking_and_normalization():
    p = np.array([[0.0, 0.1, 0.9], [0.2, 0.4, 0.4]])
    result = dev.temperature_scale(p, 1.5)
    assert np.isfinite(result).all()
    np.testing.assert_allclose(result.sum(axis=1), 1)
    np.testing.assert_array_equal(result.argmax(axis=1), p.argmax(axis=1))
    assert result[0, 2] < 0.9


def test_sibling_pooling_preserves_strong_calls_and_requires_complete_link():
    rows = [row(a, b) for a, b in combinations("abcd", 2)]
    index = {dev.key(r): i for i, r in enumerate(rows)}
    probs = np.zeros((len(rows), len(dev.LABELS)))
    probs[:, 2] = 1
    probs[index[("a", "b")]] = [1, 0, 0, 0, 0, 0, 0, 0]
    probs[index[("a", "c")]] = [0, 0, 0.85, 0.15, 0, 0, 0, 0]
    result = dev.sibling_pool(rows, probs)
    assert result[index[("a", "c")], 2] == pytest.approx(0.925)
    assert np.array_equal(result[index[("b", "c")]], probs[index[("b", "c")]])
    # A high-support conflicting call must survive, even if relatives favor another label.
    probs[index[("a", "c")]] = [0, 0, 0.05, 0.95, 0, 0, 0, 0]
    result = dev.sibling_pool(rows, probs)
    assert np.array_equal(result[index[("a", "c")]], probs[index[("a", "c")]])


def test_source_or_model_mutation_fails_provenance(tmp_path):
    p = tmp_path / "input"
    p.write_bytes(b"original")
    bounds = {str(p): dev.digest(p)}
    dev.require_hashes(bounds)
    p.write_bytes(b"changed")
    with pytest.raises(ValueError, match="changed"):
        dev.require_hashes(bounds)


def test_confirmation_seeds_and_independence_are_disjoint():
    seen = set()
    for split, count in dev.COUNTS.items():
        seeds = {dev.SEED_BASE[split] + i for i in range(count)}
        assert not seen & seeds
        seen |= seeds
    assert dev.GATES["each_label_recall_min"] == 0.60
    assert dev.GATES["balanced_accuracy_min"] == 0.75


def test_protocol_gate_changes_fail_closed():
    p = {"schema": "task11-prospective-development-v4", "method_version": dev.VERSION, "gates": dict(dev.GATES)}
    assert dev.validate_protocol(p) == dev.GATES
    p["gates"]["balanced_accuracy_min"] = 0.50
    with pytest.raises(ValueError, match="gates"):
        dev.validate_protocol(p)


def test_confirmation_cannot_replay_development_genotypes(tmp_path):
    import argparse
    map_path, jar = tmp_path / "map", tmp_path / "jar"
    map_path.write_text("fixture")
    jar.write_text("fixture")
    protocol = tmp_path / "protocol.json"
    dev.dump(protocol, {"schema": "task11-prospective-development-v4", "method_version": dev.VERSION,
        "gates": dict(dev.GATES), "inputs": {str(p): dev.digest(p) for p in (map_path, jar)},
        "source_sha256": {}, "cohorts": [{"split": "confirmation", "index": 0}]})
    args = argparse.Namespace(protocol=protocol, split="confirmation", index=0,
        map=map_path, hapibd=jar, baseline=tmp_path / "old", candidate=None, out=tmp_path / "out")
    with pytest.raises(ValueError, match="cannot replay"):
        dev.run_cohort(args)
    assert not args.out.exists()


@pytest.mark.parametrize("label", dev.LABELS)
def test_new_confirmation_topologies_have_complete_independently_derived_truth(label):
    p = dev.focal_pedigree(label, (3, 4), True, True)
    assert len(p["pairs"]) == len(list(combinations(p["observed"], 2)))
    for a, b, relation in p["pairs"]:
        da, db = distances(p["parents"], a), distances(p["parents"], b)
        assert a not in db and b not in da
        common = set(da) & set(db)
        shortest = min(da[n] + db[n] for n in common)
        nearest = [n for n in common if da[n] + db[n] == shortest]
        depths = tuple(sorted((da[nearest[0]], db[nearest[0]])))
        left, right, count = dev.DEPTHS[relation]
        assert depths == (left, right)
        assert len(nearest) == count
    if label.startswith("avuncular"):
        assert len({relation for _, _, relation in p["pairs"]}) >= 3


def test_phase_scores_mask_real_ibd2_and_report_zero_exposure():
    maps = {"1": (np.array([1, 101]), np.array([0., 100.]))}
    both = {"1": [(1, 101, 1, 1), (1, 101, 2, 2)]}
    np.testing.assert_array_equal(dev.pair_haplotype_scores(both, maps, {}), np.zeros(6))
    one = {"1": [(1, 101, 1, 1)]}
    np.testing.assert_allclose(dev.pair_haplotype_scores(one, maps, {}), [1, 1, 100, 100, 1, 1])
    np.testing.assert_allclose(dev.pair_haplotype_scores(one, maps, {"1": [(1, 51)]}), [1, 1, 50, 50, .5, .5])
    np.testing.assert_array_equal(dev.pair_haplotype_scores(one, maps, {"1": [(1, 101)]}), np.zeros(6))


def test_phase_scores_endpoint_chromosome_flip_duplicate_and_split_invariance():
    maps = {chrom: (np.array([1, 101]), np.array([0., 100.])) for chrom in ("1", "2")}
    records = {"1": [(1, 41, 1, 1), (41, 101, 2, 1)], "2": [(1, 101, 1, 2)]}
    masks = {"1": [(1, 11)]}
    expected = dev.pair_haplotype_scores(records, maps, masks)
    # Homolog orientation can independently change in each chromosome/sample.
    flipped = {"1": [(a, b, 3-h1, h2) for a, b, h1, h2 in records["1"]],
               "2": [(a, b, h1, 3-h2) for a, b, h1, h2 in records["2"]]}
    swapped = {c: [(a, b, h2, h1) for a, b, h1, h2 in rows] for c, rows in records.items()}
    duplicated = {c: rows * 3 for c, rows in records.items()}
    split = {**records, "2": [(1, 61, 1, 2), (41, 101, 1, 2)]}
    for changed in (flipped, swapped, duplicated, split):
        np.testing.assert_allclose(dev.pair_haplotype_scores(changed, maps, masks), expected)
    assert expected[0] == pytest.approx(160 / 190)
    assert expected[1] == 1


def test_phase_scores_same_haplotype_fan_is_not_ibd2_and_bounds_fail():
    maps = {"1": (np.array([1, 101]), np.array([0., 100.]))}
    fan = {"1": [(1, 101, 1, 1), (1, 101, 1, 2)]}
    np.testing.assert_allclose(dev.pair_haplotype_scores(fan, maps, {}), [.5, 1, 100, 200, 1, 1])
    with pytest.raises(ValueError, match="map bounds"):
        dev.pair_haplotype_scores({"1": [(1, 102, 1, 1)]}, maps, {})


def test_phase_ten_cm_filter_unions_edges_and_short_edges_still_mask_ibd2():
    maps = {"1": (np.array([1, 101]), np.array([0., 100.]))}
    np.testing.assert_array_equal(dev.pair_haplotype_scores({"1": [(1, 10, 1, 1)]}, maps, {}), np.zeros(6))
    # Two overlapping pieces of one long edge must survive the fixed cutoff.
    np.testing.assert_allclose(dev.pair_haplotype_scores({"1": [(1, 7, 1, 1), (6, 12, 1, 1)]}, maps, {}), [1, 1, 11, 11, 1, 1])
    # A subthreshold opposite edge is still genuine IBD2 and removes that span.
    np.testing.assert_allclose(dev.pair_haplotype_scores({"1": [(1, 101, 1, 1), (1, 6, 2, 2)]}, maps, {}), [1, 1, 95, 95, 1, 1])
