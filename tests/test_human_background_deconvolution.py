import numpy as np

from scripts.fit_human_background_deconvolution import (
    combine_probabilities,
    remove_shortest_segments,
)
from scripts.fit_human_family_relationships import summarize_raw_lengths


def test_exact_segment_removal_updates_all_summaries():
    row = {"indID1": "a", "indID2": "b", "label": "siblings"}
    row.update(
        summarize_raw_lengths(
            [("1", 12.0), ("1", 22.0), ("2", 35.0), ("3", 80.0)],
            genome_cm=100.0,
        )
    )
    result = remove_shortest_segments([row], 2, minimum_cm=20, genome_cm=100.0)[0]
    assert result["segments_ge_20"] == 1
    assert result["total_cm_ge_20"] == 80.0
    assert result["relatedness_ge_20"] == 0.4
    assert result["ordered_length_1"] == 80.0
    assert result["ordered_length_2"] == 12.0


def test_marginal_probabilities_remain_normalized():
    stack = np.asarray(
        [
            [[0.8, 0.2], [0.3, 0.7]],
            [[0.4, 0.6], [0.9, 0.1]],
            [[0.2, 0.8], [0.6, 0.4]],
        ]
    )
    combined = combine_probabilities(stack, 1.5, 2.0, "marginal")
    np.testing.assert_allclose(combined.sum(axis=1), 1.0)


def test_map_selects_one_candidate_cleaning():
    stack = np.asarray([[[0.9, 0.1]], [[0.2, 0.8]], [[0.6, 0.4]]])
    combined = combine_probabilities(stack, 1.0, 10.0, "map")
    assert any(np.array_equal(combined, candidate) for candidate in stack)
