"""Software-only regressions for locked deadline wrapper; no recovery claims."""

import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
spec = importlib.util.spec_from_file_location(
    "deadline_sweep", HERE / "run_final_sweep.py"
)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def report(rows):
    return np.asarray(
        [
            [position + 1, position + 1, position + 1, 1, 1, 1, score]
            for position, score in rows
        ],
        dtype=float,
    )


def test_geometry_clips_inaccessible_and_edges():
    intervals, audit = m.final_geometry([19], 20, np.asarray([[0, 10], [30, 100]]), 100)
    assert intervals == [[0, 10]]
    assert audit[0]["final_center"] == 9
    assert audit[0]["center_adjusted"]
    assert m.final_geometry([99], 20, np.asarray([[0, 100]]), 100)[0] == [[79, 100]]


def test_nearest_tie_lower_and_exon_intersection():
    assert m.nearest(20, [[0, 11], [30, 40]]) == 10
    assert m.intersection([[5, 25], [35, 45]], [[0, 10], [20, 40]]).tolist() == [
        [5, 10],
        [20, 25],
        [35, 40],
    ]


def test_task8_variable_count_and_position_ties():
    result = m.task8_candidates(
        report([(4_000_000, 10), (100, 10), (1000, 11)]),
        report([(100, 10), (3_000_000, 1)]),
        np.asarray([[0, 5_000_000]]),
    )
    assert [r["raisd_position"] for r in result["candidates"]] == [1000, 4_000_000]
    assert result["threshold"] == 2
    assert result["fallback_reason"] is None
    tied = m.task8_candidates(
        report([(1000, 10), (100, 10)]),
        report([(100, 10), (3_000_000, 1)]),
        np.asarray([[0, 5_000_000]]),
    )
    assert tied["candidates"][0]["raisd_position"] == 100


@pytest.mark.parametrize(
    "reference, reason",
    [
        ([], "missing_reference_background"),
        ([(1, 100), (3_000_000, 50)], "no_above_threshold_callable_exonic_candidates"),
    ],
)
def test_task8_deterministic_one_peak_fallback(reference, reason):
    result = m.task8_candidates(
        report([(800, 4), (200, 4)]), report(reference), np.asarray([[0, 1000]])
    )
    assert len(result["candidates"]) == 1
    assert result["candidates"][0]["center"] == 200
    assert result["fallback_reason"] == reason
    assert result["calibrated_null"] is False


def test_task8_distant_exon_fallback_stays_callable():
    result = m.task8_candidates(
        report([(800, 4)]), report([]), np.asarray([[100_000, 100_100]])
    )
    assert result["candidates"][0]["center"] == 800
    assert m.final_geometry([800], 60000, [[0, 1000]], 1000, True)[0] == [[0, 1000]]


def test_task8_output_merge():
    assert m.final_geometry([40, 60], 30, [[0, 100]], 100, True)[0] == [[10, 90]]


def fixture_args(tmp_path):
    vcf, mask = tmp_path / "input.vcf", tmp_path / "mask.bed"
    rng = np.random.default_rng(2026112042)
    with vcf.open("w") as handle:
        handle.write(
            '##fileformat=VCFv4.2\n##contig=<ID=9,length=200000>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        )
        handle.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
            + "\t".join(f"s{i}" for i in range(10))
            + "\n"
        )
        for i in range(120):
            row = rng.integers(0, 2, 20)
            row[:2] = [0, 1]
            handle.write(
                f"9\t{i * 100 + 1}\t.\tA\tC\t.\tPASS\t.\tGT\t"
                + "\t".join(f"{row[j]}|{row[j + 1]}" for j in range(0, 20, 2))
                + "\n"
            )
    mask.write_text("9\t0\t200000\n")
    freeze = tmp_path / "source_freeze.json"
    m.dump(
        freeze,
        {
            "runner_sha256": m.sha256(HERE / "run_final_sweep.py"),
            "protocol_sha256": m.PROTOCOL_SHA256,
            "source_sha256": m.source_manifest(ROOT),
        },
    )
    return argparse.Namespace(
        task=6,
        vcf=vcf,
        accessible=mask,
        raisd=Path(__import__("sys").executable),
        source_root=ROOT,
        source_freeze=freeze,
        out_dir=tmp_path / "result",
        vcf_sha256=m.sha256(vcf),
        mask_sha256=m.sha256(mask),
        exons=None,
        exons_sha256=None,
        task7_reference=None,
    )


def test_real_small_task6_wrapper_and_output_validation(tmp_path):
    args = fixture_args(tmp_path)
    result = m.run(args)
    assert result["state"] == "COMPLETED"
    assert result["samples_diploid"] == 10
    assert result["scientifically_validated"] is False
    assert result["contig_length"] == 200000
    assert "deadline_submission.bed" in result["artifact_sha256"]
    with pytest.raises(ValueError, match="overwrite"):
        m.run(args)


def test_post_input_mutation_fails_closed(tmp_path, monkeypatch):
    args = fixture_args(tmp_path)
    original = m.snp_h12_peak

    def mutate(*params):
        answer = original(*params)
        args.accessible.write_text("9\t0\t199999\n")
        return answer

    monkeypatch.setattr(m, "snp_h12_peak", mutate)
    with pytest.raises(ValueError, match="changed"):
        m.run(args)
    assert not (args.out_dir / "completion.json").exists()
    assert json.loads((args.out_dir / "failure.json").read_text())["state"] == "FAILED"


def test_download_hash_mismatch_rejected(tmp_path):
    args = fixture_args(tmp_path)
    args.vcf_sha256 = "0" * 64
    with pytest.raises(ValueError, match="authenticated download"):
        m.run(args)


def test_task8_missing_reference_executes_declared_fallback(tmp_path, monkeypatch):
    args = fixture_args(tmp_path)
    args.task = 8
    args.exons = tmp_path / "exons.bed"
    args.exons.write_text("9\t5000\t7000\n")
    args.exons_sha256 = m.sha256(args.exons)

    def fake_raisd(command, cwd, **kwargs):
        (cwd / "RAiSD_Info.task8").write_text("| snps 120\n")
        (cwd / "RAiSD_Report.task8").write_text("6001 5901 6101 1 1 1 5\n")

    monkeypatch.setattr(m.subprocess, "run", fake_raisd)
    result = m.run(args)
    prediction = json.loads((args.out_dir / "raw" / "prediction.json").read_text())
    assert result["state"] == "COMPLETED"
    assert prediction["fallback_reason"] == "missing_reference_background"
    assert prediction["task7_completion_sha256"] is None
    assert len(result["final_intervals"]) == 1


def test_supplied_failed_reference_rejected(tmp_path):
    m.dump(tmp_path / "completion.json", {"state": "FAILED", "task": 7})
    with pytest.raises(ValueError, match="completed Task7"):
        m.checked_reference(tmp_path, {}, "abc", 200000)
