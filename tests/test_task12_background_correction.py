import copy

import numpy as np
import pytest

from scripts.fit_human_background_corrected import (
    correct_records,
    corrected_feature_names,
)
from scripts.fit_human_family_relationships import BIN_EDGES, summarize_raw_lengths
from scripts.task12_background_correction import (
    FULL_HISTOGRAM,
    LEGACY,
    METADATA_KEY,
    correct_for_model,
    correct_records_v2,
    correction_metadata,
    subtract_shortest_histogram,
)


GENOME_CM = 3346.298892538


def record(lengths):
    row = {"indID1": "a", "indID2": "b", "label": "full_4th_cousin"}
    row.update(summarize_raw_lengths([("1", value) for value in lengths], GENOME_CM))
    return row


def v2_model():
    return {
        METADATA_KEY: correction_metadata(
            FULL_HISTOGRAM, minimum_cm=5, genome_cm=GENOME_CM
        ),
        "feature_names": corrected_feature_names(5),
        "corrected_model_weight": 1.0,
    }


def test_full_histogram_removes_short_background_beyond_top_twenty():
    reference = [record([60.0])]
    observed = [record([60.0] + [6.0] * 100)]
    before = copy.deepcopy(observed)
    legacy = correct_records(observed, reference, 5, GENOME_CM)[0]
    repaired = correct_records_v2(
        observed, reference, minimum_cm=5, genome_cm=GENOME_CM
    )[0]
    # Reproduce the historical defect without rewriting legacy model semantics.
    assert legacy["length_bin_1"] == 81
    assert legacy["length_bin_6"] == 0
    assert repaired["length_bin_1"] == 0
    assert repaired["length_bin_6"] == 1
    assert repaired["segments_ge_5"] == 1
    assert repaired["total_cm_ge_5"] == 60
    assert observed == before
    assert reference == [record([60.0])]


def test_shortest_histogram_matches_full_length_oracle_over_many_counts():
    rng = np.random.default_rng(12002026)
    for _ in range(100):
        lengths = rng.uniform(0, 150, size=200).tolist()
        histogram = np.histogram(lengths, BIN_EDGES)[0].tolist()
        for threshold in (5, 10, 20, 40):
            eligible = sorted(value for value in lengths if value >= threshold)
            below = [value for value in lengths if value < threshold]
            removal = int(rng.integers(0, 240))
            oracle = np.histogram(below + eligible[removal:], BIN_EDGES)[0]
            assert (
                subtract_shortest_histogram(histogram, removal, minimum_cm=threshold)
                == oracle.tolist()
            )


def test_fractional_mass_saturation_and_subthreshold_bins():
    histogram = [7, 3, 2, 1, 0, 0, 0, 0, 0]
    corrected = subtract_shortest_histogram(histogram, 3.25)
    assert corrected == [7, 0, 1.75, 1, 0, 0, 0, 0, 0]
    assert sum(corrected[1:]) == sum(histogram[1:]) - 3.25
    assert subtract_shortest_histogram(histogram, 100) == [7, 0, 0, 0, 0, 0, 0, 0, 0]
    assert subtract_shortest_histogram(histogram, 0) == histogram
    assert histogram == [7, 3, 2, 1, 0, 0, 0, 0, 0]


def test_fractional_quantile_offset_matches_corrected_minimum_count():
    reference = [record([]) for _ in range(4)]
    observed = [record([6.0] * count) for count in (1, 2, 3, 4)]
    corrected = correct_records_v2(
        observed, reference, minimum_cm=5, genome_cm=GENOME_CM
    )
    assert [row["segments_ge_5"] for row in corrected] == [0, 0.25, 1.25, 2.25]
    for row in corrected:
        assert row["segments_ge_5"] == sum(row[f"length_bin_{i}"] for i in range(1, 9))


@pytest.mark.parametrize("invalid", [-1, float("nan"), float("inf"), True, "bad"])
def test_histogram_rejects_invalid_counts_and_offsets(invalid):
    with pytest.raises(ValueError, match="finite and nonnegative"):
        subtract_shortest_histogram([0] * 9, invalid)
    with pytest.raises(ValueError, match="finite and nonnegative"):
        subtract_shortest_histogram([invalid] + [0] * 8, 0)


def test_histogram_requires_complete_bins_and_aligned_threshold():
    with pytest.raises(ValueError, match="every length bin"):
        subtract_shortest_histogram([0] * 8, 0)
    with pytest.raises(ValueError, match="align"):
        subtract_shortest_histogram([0] * 9, 0, minimum_cm=7)


@pytest.mark.parametrize(
    "key,value",
    [
        ("length_bin_1", -1),
        ("total_cm_ge_5", float("nan")),
        ("segments_ge_5", 99),
        ("relatedness_ge_5", 0.5),
    ],
)
def test_v2_rejects_invalid_or_inconsistent_raw_measurements(key, value):
    observed = record([6.0, 60.0])
    observed[key] = value
    with pytest.raises(ValueError):
        correct_records_v2(
            [observed], [record([60.0])], minimum_cm=5, genome_cm=GENOME_CM
        )


def test_v2_rejects_impossible_length_total_empty_cohort_and_incomplete_features():
    observed = record([6.0])
    observed["total_cm_ge_5"] = 1
    observed["relatedness_ge_5"] = 1 / (2 * GENOME_CM)
    with pytest.raises(ValueError, match="threshold"):
        correct_records_v2([observed], [record([])], minimum_cm=5, genome_cm=GENOME_CM)
    with pytest.raises(ValueError, match="nonempty"):
        correct_records_v2([], [record([])], minimum_cm=5, genome_cm=GENOME_CM)
    with pytest.raises(ValueError, match="nonempty"):
        correct_records_v2([record([])], [], minimum_cm=5, genome_cm=GENOME_CM)
    with pytest.raises(ValueError, match="complete raw"):
        correct_records_v2([{}], [record([])], minimum_cm=5, genome_cm=GENOME_CM)


@pytest.mark.parametrize("genome", [0, -1, float("nan"), float("inf"), True])
def test_v2_rejects_invalid_genome_length(genome):
    with pytest.raises(ValueError):
        correct_records_v2([record([])], [record([])], minimum_cm=5, genome_cm=genome)


def test_old_artifacts_remain_legacy_and_new_metadata_routes_v2():
    observed, reference = [record([60.0] + [6.0] * 100)], [record([60.0])]
    expected = correct_records(observed, reference, 5, GENOME_CM)
    assert (
        correct_for_model(observed, reference, {}, minimum_cm=5, genome_cm=GENOME_CM)
        == expected
    )
    legacy_model = {
        METADATA_KEY: correction_metadata(LEGACY, minimum_cm=5, genome_cm=GENOME_CM)
    }
    assert (
        correct_for_model(
            observed, reference, legacy_model, minimum_cm=5, genome_cm=GENOME_CM
        )
        == expected
    )
    repaired = correct_for_model(
        observed, reference, v2_model(), minimum_cm=5, genome_cm=GENOME_CM
    )
    assert repaired == correct_records_v2(
        observed, reference, minimum_cm=5, genome_cm=GENOME_CM
    )
    assert repaired != expected


@pytest.mark.parametrize(
    "metadata",
    [
        None,
        {},
        {"version": FULL_HISTOGRAM},
        {"version": "unknown", "minimum_cm": 5, "genome_cm": GENOME_CM},
        {"version": FULL_HISTOGRAM, "minimum_cm": 20, "genome_cm": GENOME_CM},
        {"version": FULL_HISTOGRAM, "minimum_cm": 5, "genome_cm": GENOME_CM + 1},
    ],
)
def test_explicit_incompatible_metadata_never_falls_back(metadata):
    model = v2_model()
    model[METADATA_KEY] = metadata
    with pytest.raises(ValueError, match="metadata|version"):
        correct_for_model(
            [record([])], [record([])], model, minimum_cm=5, genome_cm=GENOME_CM
        )


@pytest.mark.parametrize(
    "change",
    [
        "ordered_feature",
        "reordered_features",
        "missing_feature",
        "blend",
        "boolean_blend",
    ],
)
def test_v2_rejects_incompatible_features_and_baseline_blends(change):
    model = v2_model()
    if change == "ordered_feature":
        model["feature_names"].append("ordered_length_1")
    elif change == "reordered_features":
        model["feature_names"].reverse()
    elif change == "missing_feature":
        model["feature_names"].pop()
    else:
        model["corrected_model_weight"] = True if change == "boolean_blend" else 0.5
    with pytest.raises(ValueError, match="contract|baseline"):
        correct_for_model(
            [record([])], [record([])], model, minimum_cm=5, genome_cm=GENOME_CM
        )


def test_v2_rejects_other_feature_thresholds():
    with pytest.raises(ValueError, match="only minimum_cm=5"):
        correct_records_v2(
            [record([])], [record([])], minimum_cm=20, genome_cm=GENOME_CM
        )
