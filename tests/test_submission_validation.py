import subprocess
import sys

import pytest

from ghist_submission.validate import all_pairs, validate_submission


def test_demography_valid_all_four(tmp_path):
    values = {
        1: "generations: 4500\npost_decline_fraction: 0.52\n",
        2: "mainland_population_size: 20000\nisland_population_size: 80000\ngenerations_since_split: 4000\ngenerations_since_migration: 1000\nmigration_rate: 0.0002\n",
        3: "generations: 5000\npost_growth_scaling: 2.0\n",
        4: "AFR_population_size: 16000\nEUR_population_size: 24000\ngenerations_since_split: 4500\nmigration_rate: 0.00015\n",
    }
    for challenge, value in values.items():
        path = tmp_path / f"{challenge}.yaml"
        path.write_text(value)
        assert validate_submission(path, challenge, "testing")["format_valid"]


@pytest.mark.parametrize("value", [".nan", ".inf", "-1", "true", '"5000"'])
def test_rejects_invalid_demographic_numbers(tmp_path, value):
    path = tmp_path / "prediction.yaml"
    path.write_text(f"generations: {value}\npost_decline_fraction: 0.52\n")
    with pytest.raises(ValueError):
        validate_submission(path, 1, "final")


def test_duplicate_yaml_keys_and_wrong_schema(tmp_path):
    path = tmp_path / "prediction.yaml"
    for value in (
        "generations: 5000\ngenerations: 1000\npost_decline_fraction: 0.5\n",
        "generations: 5000\ngrowth_magnitude: 2\n",
    ):
        path.write_text(value)
        with pytest.raises(ValueError):
            validate_submission(path, 3, "testing")


@pytest.mark.parametrize("challenge", [5, 6, 7, 8])
def test_final_sweep_rejects_testing_chromosome(tmp_path, challenge):
    path = tmp_path / "prediction.bed"
    path.write_text("chr7\t100\t200\n")
    with pytest.raises(ValueError, match="chromosome 9"):
        validate_submission(path, challenge, "final")
    path.write_text("chr9\t100\t200\n")
    assert validate_submission(path, challenge, "final")["format_valid"]


def test_pair_coverage_duplicates_labels_and_degree(tmp_path):
    path = tmp_path / "prediction.csv"
    path.write_text("indID1,indID2,meioses_count,relation\na,b,2,siblings\n")
    assert validate_submission(path, 9, "final", {("b", "a")})["format_valid"]
    with pytest.raises(ValueError, match="coverage"):
        validate_submission(path, 9, "final", all_pairs(["a", "b", "c"]))
    for row, error in [
        ("b,a,2,siblings", "duplicate"),
        ("a,c,3,siblings", "conflicts"),
        ("a,c,2,half_sibling", "unaccepted"),
        ("a,c,6,half_2nd_cousin", "unaccepted"),
        ("c,c,2,siblings", "self pair"),
    ]:
        path.write_text(
            "indID1,indID2,meioses_count,relation\na,b,2,siblings\n" + row + "\n"
        )
        with pytest.raises(ValueError, match=error):
            validate_submission(path, 9, "testing")


@pytest.mark.parametrize("challenge", [10, 11, 12])
def test_relation_only_contract(tmp_path, challenge):
    path = tmp_path / "prediction.tsv"
    path.write_text("indID1\tindID2\trelation\na\tb\tsiblings\n")
    assert validate_submission(path, challenge, "final", {("a", "b")})["format_valid"]


def test_cli_returns_failure_for_wrong_dataset(tmp_path):
    path = tmp_path / "prediction.bed"
    path.write_text("chr7\t10\t20\n")
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "scripts.validate_submission",
            str(path),
            "--challenge",
            "5",
            "--dataset",
            "final",
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 1
    assert "chromosome 9" in result.stdout
