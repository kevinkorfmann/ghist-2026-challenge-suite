import gzip
from itertools import combinations
from pathlib import Path

import numpy as np
import pytest

from scripts.reconstruct_bonsai_pedigree import (
    HaplotypeSegment,
    convert_segments,
    coverage_intervals,
    numeric_id,
    read_haplotype_segments,
    read_map,
    segment_summaries,
)


def test_numeric_id() -> None:
    assert numeric_id("id_0") == 1
    assert numeric_id("sample_41") == 42


@pytest.mark.parametrize("edges", [((1, 1), (2, 2)), ((1, 2), (2, 1))])
def test_coverage_intervals_recovers_full_sharing(edges) -> None:
    assert coverage_intervals([(10, 30, *edges[0]), (20, 40, *edges[1])]) == [
        (10, 20, False),
        (20, 30, True),
        (30, 40, False),
    ]


@pytest.mark.parametrize(
    "edges", [((1, 1), (1, 2)), ((1, 1), (2, 1)), ((1, 2), (1, 2))]
)
def test_same_haplotype_or_duplicate_edges_do_not_create_ibd2(edges) -> None:
    assert coverage_intervals([(10, 30, *edges[0]), (20, 40, *edges[1])]) == [
        (10, 40, False)
    ]


def test_all_edge_sets_require_a_two_edge_matching() -> None:
    edges = ((1, 1), (1, 2), (2, 1), (2, 2))
    for size in range(5):
        for chosen in combinations(edges, size):
            expected_full = any(
                left[0] != right[0] and left[1] != right[1]
                for left, right in combinations(chosen, 2)
            )
            expected = [(10, 20, expected_full)] if chosen else []
            assert coverage_intervals([(10, 20, *edge) for edge in chosen]) == expected


def test_nested_repeated_edges_remain_active_until_last_end() -> None:
    segments = [(10, 40, 1, 1), (10, 40, 1, 1), (15, 25, 1, 1), (20, 35, 2, 2)]
    expected = [(10, 20, False), (20, 35, True), (35, 40, False)]
    assert coverage_intervals(segments) == expected
    assert coverage_intervals(list(reversed(segments))) == expected


def test_touching_edges_do_not_create_positive_length_ibd2() -> None:
    assert coverage_intervals([(10, 20, 1, 1), (20, 30, 2, 2)]) == [(10, 30, False)]
    assert coverage_intervals([(10, 20, 1, 1), (21, 30, 2, 2)]) == [
        (10, 20, False),
        (21, 30, False),
    ]


@pytest.mark.parametrize(
    "segment",
    [
        (10, 20),
        (0, 20, 1, 1),
        (-1, 20, 1, 1),
        (10, 10, 1, 1),
        (20, 10, 1, 1),
        (float("nan"), 20, 1, 1),
        (10, float("inf"), 1, 1),
        (10.5, 20, 1, 1),
        (True, 20, 1, 1),
        (10, 20, 0, 1),
        (10, 20, 1, 3),
        (10, 20, 1.0, 2),
    ],
)
def test_coverage_rejects_invalid_or_ambiguous_segments(segment) -> None:
    with pytest.raises(ValueError):
        coverage_intervals([segment])


def test_read_plink_map(tmp_path: Path) -> None:
    path = tmp_path / "map.txt"
    path.write_text("1 . 0.0 10\n1 . 1.5 20\n")
    result = read_map(path)
    assert result["1"][0].tolist() == [10.0, 20.0]
    assert result["1"][1].tolist() == [0.0, 1.5]


def test_segment_reader_filters_to_family(tmp_path: Path) -> None:
    path = tmp_path / "segments.ibd.gz"
    with gzip.open(path, "wt") as handle:
        handle.write("id_0 1 id_1 2 1 10 20 1.0\n")
        handle.write("id_0 1 id_2 2 1 10 20 1.0\n")
    grouped, samples = read_haplotype_segments(path, {"id_0", "id_1"})
    assert set(grouped) == {("id_0", "id_1", "1")}
    assert grouped[("id_0", "id_1", "1")] == [HaplotypeSegment(10, 20, 1, 2)]
    assert samples == {"id_0", "id_1"}


def write_inputs(tmp_path: Path, rows: str, map_text: str = "1 . 0.0 10\n1 . 3.0 40\n"):
    path = tmp_path / "segments.ibd.gz"
    with gzip.open(path, "wt") as handle:
        handle.write(rows)
    map_path = tmp_path / "map.txt"
    map_path.write_text(map_text)
    return path, map_path


def test_pair_reversal_swaps_haplotypes_and_preserves_conversion(
    tmp_path: Path,
) -> None:
    path, map_path = write_inputs(
        tmp_path, "id_1 2 id_0 1 1 10 30 2\nid_0 2 id_1 1 1 20 40 2\n"
    )
    grouped, _ = read_haplotype_segments(path)
    assert grouped[("id_0", "id_1", "1")][0] == HaplotypeSegment(10, 30, 1, 2)
    converted, bio, identifiers = convert_segments(path, map_path, 0)
    assert converted == [
        [1, 2, "1", 10, 20, False, 1.0],
        [1, 2, "1", 20, 30, True, 1.0],
        [1, 2, "1", 30, 40, False, 1.0],
    ]
    reversed_path, _ = write_inputs(
        tmp_path, "id_0 1 id_1 2 1 10 30 2\nid_1 1 id_0 2 1 20 40 2\n"
    )
    assert convert_segments(reversed_path, map_path, 0) == (converted, bio, identifiers)


def test_reversed_duplicate_does_not_create_ibd2(tmp_path: Path) -> None:
    path, map_path = write_inputs(
        tmp_path, "id_1 2 id_0 1 1 10 30 2\nid_0 1 id_1 2 1 20 40 2\n"
    )
    assert convert_segments(path, map_path, 0)[0] == [[1, 2, "1", 10, 40, False, 3.0]]


def test_threshold_before_scaling_interpolation_and_zero_pairs(tmp_path: Path) -> None:
    path, map_path = write_inputs(
        tmp_path, "id_0 1 id_1 1 1 15 30 99\n", "1 . 0 10\n1 . 1 20\n1 . 3 40\n"
    )
    segments, bio, identifiers = convert_segments(
        path, map_path, 1.5, {"id_0", "id_1", "id_2"}, 2
    )
    assert segments == [[1, 2, "1", 15, 30, False, 3.0]]
    assert identifiers == {1: "id_0", 2: "id_1", 3: "id_2"}
    assert [row["genotype_id"] for row in bio] == [1, 2, 3]
    assert convert_segments(path, map_path, 2, cm_scale=2)[0] == []


def test_zero_cm_plateau_is_valid_but_not_emitted(tmp_path: Path) -> None:
    path, map_path = write_inputs(
        tmp_path, "id_0 1 id_1 1 1 10 20 0\n", "1 . 0 10\n1 . 0 20\n1 . 1 40\n"
    )
    assert convert_segments(path, map_path, 0)[0] == []


@pytest.mark.parametrize(
    "row",
    [
        "id_0 1 id_1 1 1 10 20",
        "id_0 1 id_1 1 1 nan 20 1",
        "id_0 1 id_1 1 1 10 inf 1",
        "id_0 1 id_1 1 1 0 20 1",
        "id_0 1 id_1 1 1 -1 20 1",
        "id_0 1 id_1 1 1 10 10 1",
        "id_0 1 id_1 1 1 20 10 1",
        "id_0 3 id_1 1 1 10 20 1",
        "id_0 1 id_1 0 1 10 20 1",
        "id_0 1 id_1 1 1 10 20 nan",
        "id_0 1 id_1 1 1 10 20 -1",
        "id_0 1 id_0 2 1 10 20 1",
    ],
)
def test_segment_reader_rejects_malformed_records(tmp_path: Path, row: str) -> None:
    path, _ = write_inputs(tmp_path, row + "\n")
    with pytest.raises(ValueError, match=":1:"):
        read_haplotype_segments(path)


@pytest.mark.parametrize(
    "map_text",
    [
        "",
        "1 . 0 10\n",
        "1 . 0 10\n1 . 1 10\n",
        "1 . 0 20\n1 . 1 10\n",
        "1 . 1 10\n1 . 0 20\n",
        "1 . nan 10\n1 . 1 20\n",
        "1 . 0 nan\n1 . 1 20\n",
        "1 . 0 inf\n1 . 1 20\n",
        "1 . -1 10\n1 . 1 20\n",
        "1 . 0 -1\n1 . 1 20\n",
        "1 . 0 0\n1 . 1 20\n",
        "1 . 0 10.5\n1 . 1 20\n",
        "1 . 10\n",
    ],
)
def test_map_rejects_invalid_coordinates(tmp_path: Path, map_text: str) -> None:
    _, path = write_inputs(tmp_path, "", map_text)
    with pytest.raises(ValueError):
        read_map(path)


@pytest.mark.parametrize(
    "row, message",
    [
        ("id_0 1 id_1 1 2 10 20 1\n", "missing from"),
        ("id_0 1 id_1 1 1 5 20 1\n", "exceeds"),
        ("id_0 1 id_1 1 1 10 50 1\n", "exceeds"),
    ],
)
def test_conversion_rejects_unmapped_intervals(
    tmp_path: Path, row: str, message: str
) -> None:
    path, map_path = write_inputs(tmp_path, row)
    with pytest.raises(ValueError, match=message):
        convert_segments(path, map_path, 0)


@pytest.mark.parametrize(
    "threshold, scale",
    [(-1, 1), (np.nan, 1), (np.inf, 1), (0, 0), (0, -1), (0, np.nan), (0, np.inf)],
)
def test_conversion_rejects_invalid_threshold_or_scale(
    tmp_path: Path, threshold, scale
) -> None:
    path, map_path = write_inputs(tmp_path, "")
    with pytest.raises(ValueError):
        convert_segments(path, map_path, threshold, cm_scale=scale)


def test_segment_summaries_include_zero_pairs() -> None:
    segments = [[1, 2, "24", 10, 20, False, 12.0]]
    result = segment_summaries(segments, [1, 2, 3])
    assert result[frozenset((1, 2))]["total_half"] == 12.0
    assert result[frozenset((1, 3))]["total_half"] == 0.0
