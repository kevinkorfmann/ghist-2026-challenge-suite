from pathlib import Path

from scripts.reconstruct_bonsai_pedigree import (
    coverage_intervals,
    numeric_id,
    read_haplotype_segments,
    read_map,
    segment_summaries,
)


def test_numeric_id() -> None:
    assert numeric_id("id_0") == 1
    assert numeric_id("sample_41") == 42


def test_coverage_intervals_recovers_full_sharing() -> None:
    assert coverage_intervals([(10, 30), (20, 40)]) == [
        (10, 20, False),
        (20, 30, True),
        (30, 40, False),
    ]


def test_read_plink_map(tmp_path: Path) -> None:
    path = tmp_path / "map.txt"
    path.write_text("1 . 0.0 10\n1 . 1.5 20\n")
    result = read_map(path)
    assert result["1"][0].tolist() == [10.0, 20.0]
    assert result["1"][1].tolist() == [0.0, 1.5]


def test_segment_reader_filters_to_family(tmp_path: Path) -> None:
    import gzip

    path = tmp_path / "segments.ibd.gz"
    with gzip.open(path, "wt") as handle:
        handle.write("id_0 1 id_1 2 1 10 20 1.0\n")
        handle.write("id_0 1 id_2 2 1 10 20 1.0\n")
    grouped, samples = read_haplotype_segments(path, {"id_0", "id_1"})
    assert set(grouped) == {("id_0", "id_1", "1")}
    assert samples == {"id_0", "id_1"}


def test_segment_summaries_include_zero_pairs() -> None:
    segments = [[1, 2, "24", 10, 20, False, 12.0]]
    result = segment_summaries(segments, [1, 2, 3])
    assert result[frozenset((1, 2))]["total_half"] == 12.0
    assert result[frozenset((1, 3))]["total_half"] == 0.0
