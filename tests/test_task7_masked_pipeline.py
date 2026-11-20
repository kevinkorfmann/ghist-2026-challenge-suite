"""Synthetic contract checks, not biological recovery evidence."""

import json
from pathlib import Path

import numpy as np
import pytest

from ghist_single_sweep import background
from ghist_single_sweep.scan import VariantMatrix, accessible_bases
from scripts import freeze_task7_sweep_bgs as freeze
from scripts import predict_task7_sweep_bgs as cli


def matrix(positions, length=2_000_000):
    positions = np.asarray(positions, dtype=np.int64)
    return VariantMatrix(
        "7", length, positions, np.tile([0, 0, 1, 1], (len(positions), 1)), 0
    )


def bounded_scans(monkeypatch, anchor=750_000):
    observed = []

    def scan(variants, mask, width, step):
        observed.append(variants.positions.tolist())
        assert step == width // 5
        return {
            "center": np.array([anchor, anchor + 1000]),
            "composite_score": np.array([2.0, 1.0]),
        }

    def ld(variants, candidates, width):
        return {
            "candidate_position": candidates,
            "flank_ld_ratio": np.ones(len(candidates)),
        }

    monkeypatch.setattr(background, "scan_windows", scan)
    monkeypatch.setattr(background, "flank_ld_scan", ld)
    return observed


def test_raisd_inclusive_span_conversion_and_holes():
    report = np.array([[15, 10, 20, 1, 1, 1, 2], [35, 30, 40, 1, 1, 1, 3]], dtype=float)
    # VCF10..20 is exactly BED[9,20); the first base cannot be discarded.
    assert (
        len(background.callable_raisd_windows(report[:1], np.array([[9, 20]]), 100))
        == 1
    )
    with pytest.raises(ValueError, match="No fully callable"):
        background.callable_raisd_windows(report[:1], np.array([[10, 20]]), 100)
    with pytest.raises(ValueError, match="No fully callable"):
        background.callable_raisd_windows(
            report[:1], np.array([[9, 14], [15, 20]]), 100
        )
    with pytest.raises(ValueError, match="coordinates"):
        background.callable_raisd_windows(
            report[:1] * np.nan, np.array([[0, 100]]), 100
        )


def test_gap_and_ld_masking_and_mechanical_region(monkeypatch):
    observed = bounded_scans(monkeypatch)
    variants = matrix(
        [
            100,
            1000,
            399999,
            450000,
            500000,
            600000,
            650000,
            700000,
            900000,
            1000000,
            1900000,
        ]
    )
    mask = np.array([[0, 400000], [600000, 2000000]])
    report = np.array(
        [[750001, 745001, 755001, 1, 1, 1, 10], [500001, 399999, 600001, 1, 1, 1, 100]]
    )
    result = background.predict_masked_ten(variants, mask, report)
    assert result["localization"] == {
        "region_source": "full_scan_composite_maximum",
        "region": [0, 1500000],
        "truth_used_for_inference": False,
    }
    assert result["estimates"]["largest_gap_midpoint"] == (1000 + 399999) / 2
    assert result["estimates"]["raisd"] == 750000
    assert result["raisd_windows_callable"] == 1
    assert all(
        450000 not in positions and 500000 not in positions for positions in observed
    )
    for width in (50000, 100000, 250000):
        center = result["estimates"][f"flank_ld_{width}"]
        assert (
            accessible_bases(
                np.array([center - width, center]),
                np.array([center, center + width]),
                mask,
            ).min()
            >= 0.8 * width
        )
    assert result["signal_ratio"] is None
    assert result["consensus_eligible"] is False
    assert len(result["estimates"]) == 10


def test_boundary_interval_clipped(monkeypatch):
    bounded_scans(monkeypatch, anchor=500)
    result = background.predict_masked_ten(
        matrix(np.arange(0, 1000000, 1000), 1000000),
        np.array([[0, 1000000]]),
        np.array([[501, 1, 1001, 1, 1, 1, 1]]),
    )
    assert result["center"] == 500
    assert result["interval"] == [0, 80500]


@pytest.fixture
def candidate(tmp_path, monkeypatch):
    bounded_scans(monkeypatch, anchor=500)
    vcf, mask, binary = (
        tmp_path / "input.vcf",
        tmp_path / "mask.bed",
        tmp_path / "raisd",
    )
    cli.write_masked_vcf(matrix(np.arange(0, 1000000, 1000), 1000000), vcf)
    mask.write_text("7\t0\t1000000\n")
    binary.write_bytes(b"synthetic test stub, never executed")

    def reference(command, *, cwd, **kwargs):
        Path(cwd, "RAiSD_Report.task7").write_text("501 1 1001 1 1 1 1\n")
        Path(cwd, "RAiSD_Info.task7").write_text("| snps 1000 |\n")

    monkeypatch.setattr(cli.subprocess, "run", reference)
    out = tmp_path / "output"
    result = cli.run(vcf, mask, "testing", binary, out)
    return out / "prediction.json", result


def test_cli_roundtrip_and_candidate_recompute(candidate):
    path, payload = candidate
    assert path.with_name("candidate.bed").read_text() == "7\t0\t80500\n"
    assert freeze.verify_candidate(path) == json.loads(json.dumps(payload))
    payload["center"] += 1
    path.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="recomputed"):
        freeze.verify_candidate(path)


def test_reject_stale_source_and_modified_observation(candidate):
    path, payload = candidate
    payload["source_sha256"]["src/ghist_single_sweep/background.py"] = "0" * 64
    path.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="Stale"):
        freeze.verify_candidate(path)
    payload["source_sha256"] = cli.source_manifest()
    path.write_text(json.dumps(payload))
    path.with_name("RAiSD_Report.task7").write_text("1001 1 2001 1 1 1 9\n")
    with pytest.raises(ValueError, match="artifacts"):
        freeze.verify_candidate(path)


def test_old_and_development_recovery_rejected():
    with pytest.raises(ValueError, match="Legacy"):
        freeze.validate_matched_recovery(
            {"inference_method": "median_multiscale_raisd_gap_ld", "records": []}
        )
    with pytest.raises(ValueError, match="independent confirmation"):
        freeze.validate_matched_recovery(
            {
                "schema": freeze.SCHEMA,
                "status": "development",
                "public_testing_answers_used": False,
            }
        )


def test_reference_count_mismatch(tmp_path):
    info = tmp_path / "info"
    info.write_text("| snps 49 |\n")
    with pytest.raises(ValueError, match="SNP counts"):
        cli.validate_raisd_site_count(info, 50)


def test_complete_seed_groups_required(candidate, tmp_path):
    path, prediction = candidate
    design = {
        "status": "prespecified",
        "method": background.METHOD,
        "settings": background.SETTINGS,
        "source_sha256": cli.source_manifest(),
        "simulation_context": {
            "background_selection": True,
            "slim_scaling_factor": 1,
            "demography": "explicit_stress_model",
            "dfe": "specified_DFE",
            "simulation_source_sha256": {"generator": "1" * 64},
        },
        "exposure": {
            k: prediction[k]
            for k in ("chromosome", "contig_length", "haplotypes", "callable_bases")
        },
        "software": prediction["software"],
        "calibration_seeds": [1],
        "evaluation_seeds": [1],
    }
    design["exposure"]["canonical_mask_sha256"] = prediction["artifact_sha256"][
        "accessible.bed"
    ]
    design_path = tmp_path / "design.json"
    design_path.write_text(json.dumps(design))
    recovery = {
        "schema": freeze.SCHEMA,
        "status": "complete_independent_confirmation",
        "public_testing_answers_used": False,
        "design_path": design_path.name,
        "design_sha256": cli.sha256(design_path),
        "records": [],
    }
    recovery_path = tmp_path / "recovery.json"
    recovery_path.write_text(json.dumps(recovery))
    with pytest.raises(ValueError, match="disjoint"):
        freeze.score_recovery(recovery_path, prediction)
    design["evaluation_seeds"] = [2]
    design_path.write_text(json.dumps(design))
    recovery["design_sha256"] = cli.sha256(design_path)
    recovery_path.write_text(json.dumps(recovery))
    with pytest.raises(ValueError, match="every planned seed"):
        freeze.score_recovery(recovery_path, prediction)


def test_real_small_matrix_has_finite_ten_estimator_output():
    # Actual numerical kernels; artificial genotypes provide only software coverage.
    rng = np.random.default_rng(4817)
    positions = np.arange(1000, 600000, 1000, dtype=np.int64)
    haplotypes = rng.integers(0, 2, (len(positions), 12), dtype=np.uint8)
    haplotypes[:, 0], haplotypes[:, 1] = 0, 1
    variants = VariantMatrix("7", 600000, positions, haplotypes, 0)
    result = background.predict_masked_ten(
        variants,
        np.array([[0, 600000]]),
        np.array([[300001, 275001, 325001, 1, 1, 1, 1]]),
    )
    assert len(result["estimates"]) == 10
    assert np.isfinite(list(result["estimates"].values())).all()
    assert 0 <= result["interval"][0] < result["interval"][1] <= 600000
    assert result["consensus_eligible"] is False  # no independent distant baseline


@pytest.mark.parametrize("phase_sets", [("10", "20"), ("10", "."), (".", "20")])
def test_multiple_or_mixed_phase_sets_rejected(tmp_path, phase_sets):
    path = tmp_path / "phased.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=7,length=1000>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\n"
        + "".join(
            f"7\t{10 + i}\t.\tA\tC\t.\tPASS\t.\tGT:PS\t0|1:{ps}\t0|1:1\n"
            for i, ps in enumerate(phase_sets)
        )
    )
    with pytest.raises(ValueError, match="mixed or multiple PS"):
        cli.read_phased_biallelic_vcf(path)
    path.write_text(
        path.read_text().replace(f"0|1:{phase_sets[1]}\t", f"0/0:{phase_sets[1]}\t")
    )
    assert cli.read_phased_biallelic_vcf(path).positions.tolist() == [9, 10]


def test_candidate_cannot_change_dataset_and_bed_together(candidate):
    path, payload = candidate
    payload["dataset"] = "final"
    bed = path.with_name("candidate.bed")
    bed.write_text(bed.read_text().replace("7\t", "9\t"))
    payload["artifact_sha256"]["candidate.bed"] = cli.sha256(bed)
    path.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="VCF chromosome"):
        freeze.verify_candidate(path)


@pytest.mark.parametrize("changed", ["input.vcf", "mask.bed", "raisd", "source"])
def test_input_binary_and_source_bound_before_work(
    candidate, tmp_path, monkeypatch, changed
):
    reference = cli.subprocess.run
    manifest = cli.source_manifest()

    def mutate(command, **kwargs):
        reference(command, **kwargs)
        if changed == "source":
            monkeypatch.setattr(
                cli, "source_manifest", lambda: {**manifest, "changed": "0" * 64}
            )
        else:
            with (tmp_path / changed).open("a") as handle:
                handle.write("\n")

    monkeypatch.setattr(cli.subprocess, "run", mutate)
    out = tmp_path / "changed_output"
    with pytest.raises(ValueError, match="changed during inference"):
        cli.run(
            tmp_path / "input.vcf",
            tmp_path / "mask.bed",
            "testing",
            tmp_path / "raisd",
            out,
        )
    assert not (out / "prediction.json").exists()
