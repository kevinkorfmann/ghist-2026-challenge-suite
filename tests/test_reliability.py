from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

import numpy as np
import pytest

from ghist_growth_bgs.selection import compare_growth_models
from ghist_growth_bgs.spectrum import (
    fit_growth_spectrum,
    growth_spectrum,
    sudden_growth_spectrum,
)
from ghist_secondary_contact.spectrum import fit_secondary_contact_spectrum
from ghist_secondary_contact.vcf_fast import scan_secondary_contact_vcf
from ghist_single_sweep.multiple import calibrate_maximum, detect_regions
from ghist_split_migration_bgs.spectrum import fit_isolation_with_migration_spectrum
from scripts.freeze_task8_multiple_sweeps import validated_candidate


@pytest.mark.parametrize(
    "fit,shape",
    [
        (fit_growth_spectrum, (9,)),
        (fit_secondary_contact_spectrum, (7, 5)),
        (fit_isolation_with_migration_spectrum, (7, 5)),
    ],
)
@pytest.mark.parametrize(
    "problem", ["empty", "negative", "nan", "rate", "restarts", "projection"]
)
def test_demographic_fit_rejects_invalid_data_before_optimization(fit, shape, problem):
    spectrum = np.ones(shape)
    kwargs = dict(sequence_length=10000, mutation_rate=1e-8, restarts=2, workers=1)
    if problem == "empty":
        spectrum.fill(0)
    elif problem == "negative":
        spectrum.flat[1] = -1
    elif problem == "nan":
        spectrum.flat[1] = np.nan
    elif problem == "rate":
        kwargs["mutation_rate"] = np.inf
    elif problem == "restarts":
        kwargs["restarts"] = 0
    elif problem == "projection":
        kwargs["projection"] = 100 if len(shape) == 1 else (100, 3)
    with pytest.raises(ValueError):
        fit(spectrum, **kwargs)


@pytest.mark.parametrize(
    "convention,model",
    [("exponential", growth_spectrum), ("sudden", sudden_growth_spectrum)],
)
def test_growth_fit_recovers_independent_expected_spectrum(convention, model):
    # A different parameter point than the historical rounded target.
    ratio, duration, ancestral, length, mu = 3.7, 0.31, 12345, 100000, 3.36e-8
    spectrum = np.asarray(model((ratio, duration), [12])) * (
        4 * ancestral * length * mu
    )
    fit = fit_growth_spectrum(
        spectrum,
        sequence_length=length,
        mutation_rate=mu,
        restarts=4,
        maxiter=200,
        seed=41234,
        convention=convention,
    )
    assert fit["post_growth_scaling"] == pytest.approx(ratio, rel=0.02)
    assert fit["generations"] == pytest.approx(2 * ancestral * duration, rel=0.02)


def test_chromosome_holdout_never_fits_held_out_group(monkeypatch):
    import ghist_growth_bgs.selection as selection

    spectra = np.array([[0, k, 5, 3, 2, 0] for k in (11, 22, 33)], dtype=float)
    seen = []

    def fit(data, **kwargs):
        seen.append(data.copy())
        return {
            "post_growth_scaling": 2.3,
            "scaled_duration_2N": 0.2,
            "generations": 4321.123,
        }

    monkeypatch.setattr(selection, "fit_growth_spectrum", fit)
    result = compare_growth_models(
        spectra, [100, 100, 100], ["a", "b", "c"], restarts=1
    )
    for i in range(3):
        for data in seen[2 * i : 2 * i + 2]:
            np.testing.assert_array_equal(data, spectra.sum(axis=0) - spectra[i])
    assert result["prediction"]["generations"] == 4321.123
    with pytest.raises(ValueError, match="three independent"):
        compare_growth_models(spectra, [100] * 3, ["a"] * 3)


def test_secondary_contact_is_invariant_to_vcf_sample_order(tmp_path):
    names = ["mainland_a", "mainland_b", "island_a"]
    genotypes = ["0/1", "1|1", "0/0"]
    spectra = []
    for order in ([0, 1, 2], [2, 0, 1]):
        path = tmp_path / f"{order[0]}.vcf.gz"
        with gzip.open(path, "wt") as handle:
            handle.write(
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                + "\t".join(names[i] for i in order)
                + "\n"
            )
            handle.write(
                "1\t10\t.\tA\tG\t.\tPASS\tAA=A\tGT\t"
                + "\t".join(genotypes[i] for i in order)
                + "\n"
            )
        spectrum, metadata = scan_secondary_contact_vcf(path)
        spectra.append(spectrum)
        assert spectrum.shape == (5, 3)
        assert spectrum[3, 0] == 1
        assert metadata["fully_phased_fraction"] == pytest.approx(1 / 3)
    np.testing.assert_array_equal(*spectra)


@pytest.mark.parametrize("peaks", [0, 2, 7])
def test_multiple_sweeps_has_no_preset_count_and_clips_to_chromosome(peaks):
    report = np.zeros((10, 7))
    report[:, 0] = np.arange(10) * 100 + 1
    report[:peaks, 6] = 10
    result = detect_regions(
        report,
        threshold=5,
        separation=20,
        half_width=20,
        exons=[(0, 1000)],
        chromosome_length=1000,
        max_exon_snap=0,
    )
    assert result["peak_count"] == peaks
    assert result["region_count"] == peaks
    assert all(0 <= a < b <= 1000 for a, b in result["intervals"])
    if peaks:
        assert result["peaks"][0]["center"] == 0  # VCF POS 1 is BED coordinate 0.


def test_independent_null_threshold_rejects_insufficient_evidence():
    with pytest.raises(ValueError, match="Too few"):
        calibrate_maximum([1, 2, 3], 0.05)
    assert calibrate_maximum(range(1, 20), 0.05) == 19
    with pytest.raises(ValueError):
        calibrate_maximum([np.nan] * 20)
    with pytest.raises(ValueError, match="exact candidate"):
        validated_candidate(
            {"intervals": [[1, 2]]}, {"robustness_gate": {"passed": True}}, "abc"
        )


def test_task4_export_keeps_fitted_values_instead_of_preset_targets(
    tmp_path, monkeypatch
):
    from scripts import freeze_split_migration_bgs as freeze

    directory = Path(__file__).resolve().parents[1] / "experiments/split_migration_bgs"
    args = ["freeze", "--out", str(tmp_path / "freeze.json")]
    for flag, name in (
        ("baseline", "baseline_summary"),
        ("refined", "refined_middle_p50_40"),
        ("bootstrap", "bootstrap_1mb_summary"),
        ("calibration", "selection_calibration"),
        ("near-recovery", "selection_near_observed_aggregate"),
    ):
        args += ["--" + flag, str(directory / (name + ".json"))]
    for flag in (
        "testing-vcf-sha256",
        "coding-bed-sha256",
        "recombination-map-sha256",
        "final-sealed-sha256",
    ):
        args += ["--" + flag, "0" * 64]
    monkeypatch.setattr(sys, "argv", args)
    with pytest.raises(ValueError, match="physical-unit"):
        freeze.main()
    # Isolate continuous export behavior using explicitly marked test fixtures;
    # historical calibration files above must remain unmodified and rejected.
    from ghist_growth_bgs.simulation import SIMULATION_SCALING_SCHEME
    original_load = freeze.load
    def marked_fixture(path):
        result = original_load(path)
        result["simulation_scaling_scheme"] = SIMULATION_SCALING_SCHEME
        return result
    monkeypatch.setattr(freeze, "load", marked_fixture)
    freeze.main()
    result = json.loads((tmp_path / "freeze.json").read_text())
    raw = json.loads((directory / "baseline_summary.json").read_text())[
        "projection_median"
    ]
    assert result["testing_prediction"] == {k: raw[k] for k in freeze.PARAMETERS}
    assert result["testing_prediction"]["AFR_population_size"] != 16000


@pytest.mark.parametrize(
    "position,genotype,error",
    [
        (0, "0|1", "1-based"),
        (10, "0|.", "Missing genotype"),
        (10, "0|2", "Missing genotype|biallelic diploid"),
    ],
)
def test_growth_extractors_reject_biased_or_invalid_fixed_sample_spectra(
    tmp_path, position, genotype, error
):
    from ghist_growth_bgs.audit import audit_growth_resources
    from ghist_growth_bgs.blocks import extract_block_spectra

    vcf, bed, recombination = [
        tmp_path / name for name in ("input.vcf", "exons.bed", "map.tsv")
    ]
    vcf.write_text(
        '##fileformat=VCFv4.2\n##contig=<ID=2L,length=1000>\n##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral">\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n'
        + f"2L\t{position}\t.\tA\tG\t.\tPASS\tAA=A\tGT\t{genotype}\n"
    )
    bed.write_text("2L\t0\t20\n")
    recombination.write_text("Chromosome\tPosition(bp)\tRate(cM/Mb)\n2L\t0\t0\n")
    with pytest.raises(ValueError, match=error):
        audit_growth_resources(vcf, bed, recombination)
    with pytest.raises(ValueError, match=error):
        extract_block_spectra(vcf, bed)


def test_family_coverage_uses_membership_not_just_pair_counts():
    from ghist_submission.validate import within_family_pairs

    assert within_family_pairs([["b", "a"], ["c"]], ["a", "b", "c"]) == {("a", "b")}
    with pytest.raises(ValueError, match="partition"):
        within_family_pairs([["a", "b"], ["b"]], ["a", "b", "c"])
    with pytest.raises(ValueError, match="partition"):
        within_family_pairs([["a", "b"], ["x"]], ["a", "b", "c"])


@pytest.mark.parametrize(
    "position,genotype,message",
    [(0, "0|1", "1-based"), (10, "0|.", "Missing genotype")],
)
def test_bottleneck_blocks_reject_negative_index_and_incomplete_first_site(
    tmp_path, position, genotype, message
):
    from ghist_bottleneck.sfs_baseline import complete_block_spectra_from_vcf

    path = tmp_path / "bottleneck.vcf"
    path.write_text(
        '##fileformat=VCFv4.2\n##contig=<ID=1,length=1000>\n##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral">\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n'
        + f"1\t{position}\t.\tA\tG\t.\tPASS\tAA=A\tGT\t{genotype}\n"
    )
    with pytest.raises(ValueError, match=message):
        complete_block_spectra_from_vcf(path, sequence_length=1000, block_length=100)


def test_bottleneck_retains_small_valid_fraction_without_rounding_to_zero():
    from ghist_bottleneck.infer import prediction_from_latent

    prediction = prediction_from_latent([np.log(1234.567), -20])
    assert prediction["generations"] == pytest.approx(1234.567)
    assert 0 < prediction["post_decline_fraction"] < 1e-6
    with pytest.raises(ValueError):
        prediction_from_latent([1000, 0])
    with pytest.raises(ValueError):
        prediction_from_latent([2, np.nan])


def test_bottleneck_fit_survives_one_failed_restart(monkeypatch):
    import moments
    from ghist_bottleneck.sfs_baseline import fit_two_epoch_spectrum

    calls = []

    def optimizer(*args, **kwargs):
        calls.append(1)
        if len(calls) == 1:
            raise ArithmeticError("one numerical failure")
        return np.array([0.4, 0.2])

    monkeypatch.setattr(moments.Inference, "optimize_log_lbfgsb", optimizer)
    spectrum = np.asarray(moments.Demographics1D.two_epoch([0.4, 0.2], [10])) * 200
    fit = fit_two_epoch_spectrum(spectrum, restarts=2)
    assert fit["successful_restarts"] == 1
    assert len(fit["failed_restarts"]) == 1
    assert fit["post_decline_fraction"] == pytest.approx(0.4)
