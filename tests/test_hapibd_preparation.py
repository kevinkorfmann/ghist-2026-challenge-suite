from __future__ import annotations

import gzip
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pytest

from scripts.prepare_hapibd_vcf import prepare_vcf


class HapIbdPreparationTests(unittest.TestCase):
    def test_imputes_missing_alleles_with_major_allele(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "input.vcf.gz"
            output = root / "output.vcf.gz"
            with gzip.open(source, "wt") as handle:
                handle.write("##fileformat=VCFv4.2\n")
                handle.write(
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\n"
                )
                handle.write("1\t1\t.\tA\tG\t.\tPASS\t.\tGT\t1|1\t.|0\n")
                handle.write("1\t2\t.\tA\tG,T\t.\tPASS\t.\tGT\t1|2\t0|0\n")
            subprocess.run(
                [
                    sys.executable,
                    "scripts/prepare_hapibd_vcf.py",
                    "--input",
                    str(source),
                    "--output",
                    str(output),
                    "--missing-policy",
                    "major",
                ],
                check=True,
            )
            with gzip.open(output, "rt") as handle:
                records = [line for line in handle if not line.startswith("#")]
                record = records[0]
        self.assertEqual(len(records), 1)
        self.assertTrue(record.rstrip().endswith("1|1\t1|0"))


HEADER = (
    "##fileformat=VCFv4.2\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\n"
)


def write_input(path, records):
    text = HEADER + "".join(
        f"1\t{index}\t.\tA\tG\t.\tPASS\t.\t{record}\n"
        for index, record in enumerate(records, 1)
    )
    if path.suffix == ".gz":
        with gzip.open(path, "wt") as handle:
            handle.write(text)
    else:
        path.write_text(text)
    return text


@pytest.mark.parametrize("suffix", [".vcf", ".vcf.gz"])
def test_preserves_genuine_phase_and_other_format_fields(tmp_path, suffix):
    source, output = tmp_path / f"source{suffix}", tmp_path / f"output{suffix}"
    original = write_input(source, ["GT:DP\t1|0:12\t0|1:8"])
    report = prepare_vcf(source, output)
    if suffix == ".vcf.gz":
        with gzip.open(output, "rt") as handle:
            result = handle.read()
    else:
        result = output.read_text()
    assert result == original
    assert report["imputed_alleles"] == 0


@pytest.mark.parametrize("policy", ["error", "major"])
def test_unphased_heterozygotes_never_become_phased(tmp_path, policy):
    source, output = tmp_path / "source.vcf", tmp_path / "output.vcf"
    write_input(source, ["GT\t0/1\t0|0"])
    with pytest.raises(ValueError, match="unphased heterozygote"):
        prepare_vcf(source, output, missing_policy=policy)
    assert not output.exists()


def test_late_missing_call_preserves_previous_output_and_cleans_temporary(tmp_path):
    source, output = tmp_path / "source.vcf.gz", tmp_path / "output.vcf.gz"
    write_input(source, ["GT\t0|1\t1|0", "GT\t0|.\t1|1"])
    output.write_bytes(b"previous valid artifact")
    with pytest.raises(ValueError, match="missing allele"):
        prepare_vcf(source, output)
    assert output.read_bytes() == b"previous valid artifact"
    assert sorted(p.name for p in tmp_path.iterdir()) == [
        "output.vcf.gz",
        "source.vcf.gz",
    ]


def test_unphased_homozygote_normalization_does_not_change_alleles(tmp_path):
    source, output = tmp_path / "source.vcf", tmp_path / "output.vcf"
    write_input(source, ["GT:DP\t0/0:12\t1/1:9"])
    report = prepare_vcf(source, output)
    assert output.read_text().endswith("GT:DP\t0|0:12\t1|1:9\n")
    assert report["unphased_homozygotes_normalized"] == 2


def test_explicit_imputation_counts_all_missing_site_and_replacements(tmp_path):
    source, output = tmp_path / "source.vcf", tmp_path / "output.vcf"
    write_input(source, ["GT\t.|.\t.|.", "GT\t1|1\t.|0"])
    report = prepare_vcf(source, output, missing_policy="major")
    assert report["all_missing_records_skipped"] == 1
    assert report["imputed_alleles"] == 1
    assert report["records_written"] == 1
    assert output.read_text().endswith("GT\t1|1\t1|0\n")


def test_imputation_cannot_assign_phase_to_partial_unphased_call(tmp_path):
    source, output = tmp_path / "source.vcf", tmp_path / "output.vcf"
    write_input(source, ["GT\t1|1\t./0"])
    with pytest.raises(ValueError, match="imputation cannot establish phase"):
        prepare_vcf(source, output, missing_policy="major")


def test_explicit_drop_site_policy_preserves_complete_calls_without_imputation(
    tmp_path,
):
    source, output = tmp_path / "source.vcf", tmp_path / "output.vcf"
    write_input(source, ["GT\t0|1\t1|0", "GT\t.|1\t0|0", "GT\t.|.\t.|."])
    report = prepare_vcf(source, output, missing_policy="drop-site")
    assert report["records_written"] == 1
    assert report["missing_records_dropped"] == 2
    assert report["all_missing_records_skipped"] == 1
    assert report["imputed_alleles"] == 0
    assert output.read_text() == HEADER + "1\t1\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\t1|0\n"


@pytest.mark.parametrize(
    "record, message",
    [
        ("GT\t0\t0|1", "diploid"),
        ("GT\t0|2\t0|1", "biallelic"),
        ("DP\t10\t20", "GT FORMAT"),
        ("GT\t0|1", "field count"),
    ],
)
def test_rejects_unsupported_input_contracts(tmp_path, record, message):
    source, output = tmp_path / "source.vcf", tmp_path / "output.vcf"
    write_input(source, [record])
    with pytest.raises(ValueError, match=message):
        prepare_vcf(source, output)


def test_header_only_output_is_not_published(tmp_path):
    source, output = tmp_path / "source.vcf", tmp_path / "output.vcf"
    source.write_text(HEADER)
    with pytest.raises(ValueError, match="No usable genotype"):
        prepare_vcf(source, output)
    assert not output.exists()


@pytest.mark.parametrize(
    "second",
    ["GT:PS\t1|0:2\t0|1:7", "GT:PS\t1|0:.\t0|1:7", "GT\t1|0\t0|1"],
)
@pytest.mark.parametrize("reverse", [False, True])
def test_multiple_or_mixed_phase_sets_preserve_previous_output(
    tmp_path, second, reverse
):
    source, output = tmp_path / "source.vcf", tmp_path / "output.vcf"
    records = ["GT:PS\t0|1:1\t1|0:7", second]
    write_input(source, list(reversed(records)) if reverse else records)
    output.write_text("previous valid output")
    with pytest.raises(ValueError, match="mixed or multiple phase sets"):
        prepare_vcf(source, output)
    assert output.read_text() == "previous valid output"
    assert sorted(path.name for path in tmp_path.iterdir()) == [
        "output.vcf",
        "source.vcf",
    ]


def test_one_phase_set_per_sample_and_chromosome_preserves_calls(tmp_path):
    source, output = tmp_path / "source.vcf", tmp_path / "output.vcf"
    text = HEADER + (
        "1\t1\t.\tA\tG\t.\tPASS\t.\tGT:PS\t0|1:1\t1|0:7\n"
        "1\t2\t.\tA\tG\t.\tPASS\t.\tGT:PS\t1|0:1\t0|1:7\n"
        "1\t3\t.\tA\tG\t.\tPASS\t.\tGT:PS\t0|0:.\t1|1\n"
        "2\t1\t.\tA\tG\t.\tPASS\t.\tGT:PS\t0|1:3\t1|0:8\n"
    )
    source.write_text(text)
    report = prepare_vcf(source, output)
    assert output.read_text() == text
    assert report["records_written"] == 4


def test_implicit_chromosome_wide_phase_accepts_missing_ps_and_homozygotes(tmp_path):
    source, output = tmp_path / "source.vcf", tmp_path / "output.vcf"
    original = write_input(
        source,
        [
            "GT\t0|1\t1|0",
            "GT:PS\t1|0:.\t0|1",
            "GT:PS\t0|0:2\t1|1:9",
        ],
    )
    assert prepare_vcf(source, output)["records_written"] == 3
    assert output.read_text() == original


def test_dropped_missing_site_does_not_connect_phase_sets(tmp_path):
    source, output = tmp_path / "source.vcf", tmp_path / "output.vcf"
    write_input(
        source,
        [
            "GT:PS\t0|1:1\t1|0:7",
            "GT:PS\t1|0:2\t.|.:8",
            "GT:PS\t1|0:1\t0|1:7",
        ],
    )
    report = prepare_vcf(source, output, missing_policy="drop-site")
    assert report["records_written"] == 2
    assert report["missing_records_dropped"] == 1


def test_major_imputation_cannot_bridge_independent_phase_sets(tmp_path):
    source, output = tmp_path / "source.vcf", tmp_path / "output.vcf"
    write_input(source, ["GT:PS\t0|1:1\t1|0:7", "GT:PS\t.|0:2\t1|1:7"])
    with pytest.raises(ValueError, match="mixed or multiple phase sets"):
        prepare_vcf(source, output, missing_policy="major")
    assert not output.exists()


@pytest.mark.parametrize("ps", ["-1", "1.5", "nan", "2147483648"])
def test_invalid_informative_phase_set_identifier_is_rejected(tmp_path, ps):
    source, output = tmp_path / "source.vcf", tmp_path / "output.vcf"
    write_input(source, [f"GT:PS\t0|1:{ps}\t1|0:7"])
    with pytest.raises(ValueError, match="PS must be"):
        prepare_vcf(source, output)
    assert not output.exists()


if __name__ == "__main__":
    unittest.main()
