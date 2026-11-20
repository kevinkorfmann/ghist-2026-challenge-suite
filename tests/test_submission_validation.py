import json
import subprocess
import sys

import pytest

from ghist_submission.validate import (
    all_pairs,
    validate_submission,
    within_family_pairs,
)


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


def test_large_finite_yaml_integer_does_not_overflow(tmp_path):
    path = tmp_path / "prediction.yaml"
    value = 10**400
    path.write_text(f"generations: {value}\npost_decline_fraction: 0.52\n")
    assert validate_submission(path, 1, "testing")["parameters"]["generations"] == value


@pytest.mark.parametrize(
    "content",
    [
        'indID1,indID2,meioses_count\na,b,"2\n',
        'indID1,indID2,relation\na,b,"siblings"junk\n',
        'indID1,indID2,"relation\na,b,siblings\n',
    ],
)
def test_rejects_malformed_csv_quoting(tmp_path, content):
    path = tmp_path / "prediction.csv"
    path.write_text(content)
    with pytest.raises(ValueError, match="CSV"):
        validate_submission(path, 9, "testing")


@pytest.mark.parametrize(
    "expected",
    [
        [("a",)],
        [("a", "b", "c")],
        ["ab"],
        [(None, "b")],
        [("a", " ")],
        [("a", "a")],
        42,
    ],
)
def test_rejects_invalid_expected_pair_domain(tmp_path, expected):
    path = tmp_path / "prediction.csv"
    path.write_text("indID1,indID2,relation\na,b,siblings\n")
    with pytest.raises(ValueError, match="sample ID"):
        validate_submission(path, 9, "testing", expected)


def test_pair_coverage_normalizes_ids_consistently(tmp_path):
    path = tmp_path / "prediction.csv"
    path.write_text("indID1,indID2,relation\n a , b ,siblings\n")
    expected = all_pairs([" a ", "b"])
    assert expected == {("a", "b")}
    assert validate_submission(path, 9, "testing", [(" b ", "a")])[
        "input_pair_coverage_verified"
    ]


@pytest.mark.parametrize("samples", [["a", " a "], ["a", None], ["a", " "], "ab"])
def test_rejects_ambiguous_or_invalid_cohort_ids(samples):
    with pytest.raises(ValueError, match="sample ID"):
        all_pairs(samples)


def test_family_partition_accepts_single_pass_iterables():
    families = (iter(family) for family in [["b", "a"], ["c"]])
    assert within_family_pairs(families, ["a", "b", "c"]) == {("a", "b")}


def test_vcf_cohort_ids_remain_in_header_order(tmp_path):
    from ghist_submission.validate import vcf_samples

    path = tmp_path / "cohort.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tb\ta\tc\n"
    )
    assert vcf_samples(path) == ["b", "a", "c"]


@pytest.mark.parametrize(
    "families,samples",
    [
        ([["a", "b"]], ["a", "b", "b"]),
        (["ab"], ["a", "b"]),
        ([["a", "b"], []], ["a", "b"]),
        ([["a", "b"], [" b "]], ["a", "b"]),
        ([["a", "b"]], ["a", "b", None]),
    ],
)
def test_rejects_malformed_family_partitions(families, samples):
    with pytest.raises(ValueError):
        within_family_pairs(families, samples)


@pytest.mark.parametrize("challenge", [1, 5, 8])
def test_rejects_inapplicable_pair_coverage(tmp_path, challenge):
    with pytest.raises(ValueError, match="Pair coverage applies only"):
        validate_submission(tmp_path / "unused", challenge, "testing", {("a", "b")})


@pytest.mark.parametrize(
    "content,error",
    [
        ("indID1,indID2,indID1\nx,b,a\n", "unique"),
        ("indID1,indID2\na,b,ignored\n", "extra fields"),
        ("indID1,indID2\na\n", "missing"),
        ("indID1,indID2\na,\n", "nonempty"),
        ("indID1,indID2\na,a\n", "self pair"),
        ("indID1,other\na,b\n", "headers"),
        ('indID1,indID2\na,"b', "unexpected end"),
    ],
)
def test_cli_rejects_malformed_expected_pairs_with_json_error(
    tmp_path, monkeypatch, capsys, content, error
):
    from scripts.validate_submission import main

    prediction = tmp_path / "prediction.csv"
    prediction.write_text("indID1,indID2,relation\na,b,siblings\n")
    pairs = tmp_path / "pairs.csv"
    pairs.write_text(content)
    report = tmp_path / "reports" / "validation.json"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "validate_submission",
            str(prediction),
            "--challenge",
            "9",
            "--dataset",
            "testing",
            "--pairs",
            str(pairs),
            "--report",
            str(report),
        ],
    )
    with pytest.raises(SystemExit) as exit_info:
        main()
    assert exit_info.value.code == 1
    result = json.loads(capsys.readouterr().out)
    assert result["format_valid"] is False
    assert error in result["error"]
    assert json.loads(report.read_text()) == result


def test_cli_preserves_pair_lists_with_extra_named_columns_and_reversed_duplicates(
    tmp_path, monkeypatch, capsys
):
    from scripts.validate_submission import main

    prediction = tmp_path / "prediction.csv"
    prediction.write_text('indID1,indID2,relation\n"a","b","siblings"\n')
    pairs = tmp_path / "pairs.csv"
    pairs.write_text("indID1,indID2,note\nb,a,first\na,b,\n")
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "validate_submission",
            str(prediction),
            "--challenge",
            "9",
            "--dataset",
            "testing",
            "--pairs",
            str(pairs),
        ],
    )
    main()
    result = json.loads(capsys.readouterr().out)
    assert result["format_valid"] and result["input_pair_coverage_verified"]
