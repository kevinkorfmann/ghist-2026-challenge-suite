"""Observed-interval and invariance checks for the fixed context diagnostic."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

SPEC = importlib.util.spec_from_file_location("context", Path(__file__).with_name("diagnose_family_context.py"))
context = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(context)


def fixture():
    maps = {"1": (np.array([0, 100]), np.array([0., 100.]))}
    rows = [{"indID1": a, "indID2": b, "family": 0} for a, b in (("a", "b"), ("a", "c"), ("b", "c"))]
    probabilities = np.array([[0., 1.], [0., 1.], [.99, .01]])
    records = {("a", "b", "1"): [(0, 80, 1, 1)], ("a", "c", "1"): [(20, 90, 2, 1)],
               ("b", "c", "1"): [(50, 60, 1, 1)]}
    masks = {t: {("a", "1"): [(30, 40)]} for t in (2., 5., 10.)}
    return rows, probabilities, records, masks, maps


def test_interval_union_intersection_subtraction():
    assert context.merge([(0, 10), (5, 20), (20, 30), (0, 10)]) == [(0, 30)]
    assert context.intersect([(0, 10), (20, 30)], [(5, 25)]) == [(5, 10), (20, 25)]
    assert context.subtract([(0, 30)], [(5, 10), (8, 20)]) == [(0, 5), (20, 30)]
    with pytest.raises(ValueError, match="Nonpositive"):
        context.merge([(1, 1)])


def test_011_requires_siblings_ibd0_and_masks_all_three():
    _, _, records, masks, maps = fixture()
    all_ibd, eligible = context.interval_index(records, maps)
    values = context.side_context("a", {"b", "c"}, all_ibd, eligible, masks[2.], maps)
    np.testing.assert_allclose(values, [.8, .65, .4, .5, .8, .4])
    masks[2.][("c", "1")] = [(20, 30)]
    values = context.side_context("a", {"b", "c"}, all_ibd, eligible, masks[2.], maps)
    np.testing.assert_allclose(values[2:], [.3, 3/7, .7, .3])


def test_short_sibling_sharing_still_excludes_ibd0():
    _, _, records, masks, maps = fixture()
    records[("b", "c", "1")] = [(50, 55, 1, 1)]
    all_ibd, eligible = context.interval_index(records, maps)
    assert eligible[(frozenset(("b", "c")), "1")] == []
    values = context.side_context("a", {"b", "c"}, all_ibd, eligible, masks[2.], maps)
    assert values[2] == pytest.approx(.45)


def test_ambiguous_sibling_component_becomes_all_singletons():
    rows, probabilities, *_ = fixture()
    probabilities[0] = [.99, .01]
    groups, audit = context.predicted_groups(rows, probabilities, 0)
    assert groups == {s: frozenset((s,)) for s in "abc"}
    assert audit["ambiguous_components"] == 1
    with pytest.raises(ValueError, match="Incomplete"):
        context.predicted_groups(rows[:-1], probabilities[:-1], 0)


def test_endpoint_order_sample_ids_truth_fields_and_chromosome_flips_invariant():
    rows, probabilities, records, masks, maps = fixture()
    expected, _ = context.context_features(rows, probabilities, 0, records, masks, maps)
    rename = {"a": "z", "b": "x", "c": "y"}
    changed_rows = [{"indID1": rename[r["indID2"]], "indID2": rename[r["indID1"]], "family": 18,
                     "relation": "irrelevant", "pedigree": "not read"} for r in rows[::-1]]
    changed_records = {(rename[b], rename[a], chrom): [(s, e, 3-hb, 3-ha) for s, e, ha, hb in spans]
                       for (a, b, chrom), spans in records.items()}
    changed_masks = {t: {(rename[s], chrom): spans for (s, chrom), spans in values.items()} for t, values in masks.items()}
    actual, _ = context.context_features(changed_rows, probabilities[::-1], 0, changed_records, changed_masks, maps)
    np.testing.assert_allclose(actual[::-1], expected, rtol=0., atol=0.)


def test_duplicate_split_and_overlapping_edges_do_not_change_features():
    rows, probabilities, records, masks, maps = fixture()
    expected, _ = context.context_features(rows, probabilities, 0, records, masks, maps)
    records[("a", "b", "1")] = [(0, 50, 1, 1), (45, 80, 1, 1), (0, 50, 1, 1)]
    actual, _ = context.context_features(rows, probabilities, 0, records, masks, maps)
    np.testing.assert_allclose(actual, expected, rtol=0., atol=0.)


def test_zero_exposure_and_singleton_fallback_are_explicit():
    rows, probabilities, records, masks, maps = fixture()
    singleton, _ = context.context_features(rows, np.tile([0., 1.], (3, 1)), 0, records, masks, maps)
    assert np.all(singleton[:, [i*12+j for i in range(3) for j in range(4, 12)]] == 0.)
    masks = {t: {(s, "1"): [(0, 100)] for s in "abc"} for t in masks}
    empty, audit = context.context_features(rows, probabilities, 0, records, masks, maps)
    assert np.all(empty == 0.)
    assert audit["zero_ibd0_exposure_endpoints"] == [6, 6, 6]


def test_map_bounds_fail_closed():
    _, _, records, _, maps = fixture()
    records[("a", "b", "1")] = [(-1, 80, 1, 1)]
    with pytest.raises(ValueError, match="exceeds map"):
        context.interval_index(records, maps)
