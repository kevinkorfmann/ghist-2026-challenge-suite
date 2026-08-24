from pathlib import Path

from scripts.reconstruct_bonsai_pedigree import coverage_intervals, numeric_id, read_map


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
