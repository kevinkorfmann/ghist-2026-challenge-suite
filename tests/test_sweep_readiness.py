import hashlib
import json

import numpy as np
import pytest

from ghist_single_sweep import recovery
from ghist_single_sweep.export import write_biallelic_vcf
from ghist_single_sweep.scan import VariantMatrix, scan_windows
from ghist_single_sweep.scan import read_phased_biallelic_vcf
from scripts.predict_single_sweep import run
from ghist_single_sweep.soft import snp_h12_peak
from scripts.summarize_sweep_readiness import binomial_interval, summarize
from scripts.summarize_soft_sweep_observed_transfer import summarize_transfer
from scripts.summarize_soft_sweep_development import summarize as summarize_development
from ghist_single_sweep.peaks import highest_peak_position


def test_uncallable_windows_cannot_change_normalization():
    positions = np.asarray([5, 10, 15, 105, 120, 205, 210, 220, 225, 305, 320, 330])
    rows = np.asarray([[0, 0, 1, 1], [0, 1, 0, 1], [1, 0, 0, 0]], dtype=np.uint8)
    haplotypes = rows[np.arange(len(positions)) % 3]
    short = VariantMatrix("7", 400, positions, haplotypes, 0)
    masked_extension = VariantMatrix("7", 800, positions, haplotypes, 0)
    baseline = scan_windows(short, np.asarray([[0, 400]]), 100, 100)
    extended = scan_windows(masked_extension, np.asarray([[0, 400]]), 100, 100)
    for name in baseline:
        if name.startswith("z_") or name.endswith("composite_score"):
            np.testing.assert_allclose(
                baseline[name], extended[name][:4], equal_nan=True
            )


def test_localization_cannot_select_inaccessible_gap(monkeypatch):
    monkeypatch.setattr(
        recovery,
        "flank_ld_scan",
        lambda variants, candidates, flank_width: {
            "candidate_position": candidates,
            "flank_ld_ratio": np.ones(len(candidates)),
        },
    )
    positions = np.asarray([10, 20, 30, 100_000, 500_000, 500_100, 900_000])
    variants = VariantMatrix(
        "7", 1_000_000, positions, np.zeros((len(positions), 4)), 0
    )
    result = recovery.localize_from_scan(
        variants,
        {"composite_score": np.asarray([1.0]), "center": np.asarray([500_000])},
        accessibility=np.asarray([[0, 200_000], [500_000, 1_000_000]]),
    )
    assert result["largest_gap_midpoint"] == 700_050
    assert result["largest_gap_length"] == 399_900


def test_production_reader_prediction_bed_and_validator(tmp_path):
    vcf = tmp_path / "test.vcf"
    header = (
        "##fileformat=VCFv4.2\n##contig=<ID=7,length=400000>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\tc\n"
    )
    genotypes = ["0|0\t0|1\t1|1", "0|0\t1|0\t1|1", "0|1\t0|0\t1|1"]
    vcf.write_text(
        header
        + "".join(
            f"7\t{position}\t.\tA\tC\t.\tPASS\t.\tGT\t{genotypes[i % 3]}\n"
            for i, position in enumerate(range(100, 400_000, 400))
        )
    )
    bed = tmp_path / "callable.bed"
    bed.write_text("7\t0\t400000\n")
    out = tmp_path / "candidate.json"
    candidate = tmp_path / "candidate.bed"
    result = run(vcf, bed, "testing", candidate, out)
    assert result["format_validation"]["format_valid"]
    assert result["haplotypes"] == 6
    assert result["localization"]["truth_used_for_inference"] is False
    assert json.loads(out.read_text())["interval"] == result["interval"]
    with pytest.raises(ValueError, match="overwrite"):
        run(vcf, bed, "testing", candidate, out)


def test_unique_haplotypes_are_not_an_informative_soft_peak():
    variants = VariantMatrix(
        "7", 100, np.asarray([10, 20, 30, 40]), np.eye(4, dtype=np.uint8), 0
    )
    result = snp_h12_peak(variants, np.asarray([[0, 100]]), 4)
    assert not result["informative"]
    assert result["h12"] == pytest.approx(0.375)


def test_zero_events_does_not_certify_five_percent_error():
    interval = binomial_interval(0, 8)
    assert interval[0] == 0
    assert interval[1] == pytest.approx(0.36941664755)


def test_summary_retains_failed_and_unrun_cases_in_planned_denominator(tmp_path):
    design = {
        "cases": [
            {
                "case": i,
                "seed": 10 + i,
                "history": "constant",
                "selection_coefficient": 0.1,
            }
            for i in range(2)
        ]
    }
    design_path = tmp_path / "design.json"
    design_path.write_text(json.dumps(design))
    folder = tmp_path / "case_000"
    folder.mkdir()
    status = {
        "case": 0,
        "seed": 10,
        "state": "FAILED",
        "design_sha256": hashlib.sha256(design_path.read_bytes()).hexdigest(),
    }
    (folder / "status.json").write_text(json.dumps(status))
    result = summarize(design_path, tmp_path)
    assert result["planned_denominator"]["failed"] == 1
    assert result["planned_denominator"]["not_run"] == 1
    assert result["planned_denominator"]["coverage_success_yield_lower_bound"] == 0
    assert result["completed_only_metrics"] == {"completed": 0}
    assert not result["all_planned_complete"]
    status["seed"] = 999
    (folder / "status.json").write_text(json.dumps(status))
    with pytest.raises(ValueError, match="frozen design"):
        summarize(design_path, tmp_path)


def test_numeric_slim_alleles_preserve_genotypes_and_exact_contig_length(tmp_path):
    import tskit

    tables = tskit.TableCollection(100)
    for _ in range(2):
        individual = tables.individuals.add_row()
        for _ in range(2):
            tables.nodes.add_row(
                flags=tskit.NODE_IS_SAMPLE, time=0, individual=individual
            )
    root = tables.nodes.add_row(time=1)
    for node in range(4):
        tables.edges.add_row(0, 100, root, node)
    for position in (0, 99):
        site = tables.sites.add_row(position=position, ancestral_state="")
        tables.mutations.add_row(site=site, node=0, derived_state="12345,67890")
    tables.sort()
    ts = tables.tree_sequence()
    path = tmp_path / "canonical.vcf"
    metadata = write_biallelic_vcf(ts, path)
    variants = read_phased_biallelic_vcf(path)
    assert metadata["accepted_sites_with_non_nucleotide_source_labels"] == 2
    assert variants.contig_length == 100
    assert variants.positions.tolist() == [0, 99]
    np.testing.assert_array_equal(variants.haplotypes, ts.genotype_matrix())
    text = path.read_text().replace("0|0", ".|0", 1)
    path.write_text(text)
    with pytest.raises(ValueError, match="callable-exposure"):
        read_phased_biallelic_vcf(path)


def test_soft_transfer_matches_scan_scale_and_handles_empty_groups():
    record = {
        "seed": 75_226_024,
        "scan_window_size": 100_000,
        "localization": {
            "truth_used_for_inference": False,
            "region_source": "full_scan_composite_maximum",
        },
        "truth": {"position": 100, "minimum_end_frequency": 0.95},
        "estimates": {"soft_composite": 100, "h12": 100, "flank_ld_100kb": 100},
        "diagnostics": {"peak_soft_composite_score": 2},
    }
    recovered = {
        "records": [record],
        "selected_method": "median_soft_linkage",
        "interval_width_selection": {"selected_half_width": 75_000},
        "untouched_evaluation": {},
    }
    observed = {
        "scales": {
            str(size): {"top_windows": [{"soft_composite_score": score}]}
            for size, score in [(25_000, 100), (100_000, 2)]
        }
    }
    result = summarize_transfer(recovered, observed)
    assert result["signal_matched_evaluation"]["replicates"] == 1
    assert result["signal_matched_evaluation"]["coverage"] == 1
    observed["scales"]["100000"]["top_windows"][0]["soft_composite_score"] = 10
    result = summarize_transfer(recovered, observed)
    assert result["signal_matched_evaluation"]["replicates"] == 0
    assert result["signal_matched_evaluation"]["coverage"] is None
    assert not result["robustness_gate"]["passed"]
    record.pop("scan_window_size")
    with pytest.raises(ValueError, match="window sizes"):
        summarize_transfer(recovered, observed)


@pytest.mark.parametrize(
    "rows,match",
    [
        ([(10, "C"), (10, "G")], "sorted and unique"),
        ([(10, "<DEL>")], "single-nucleotide"),
        ([(10, "AC")], "single-nucleotide"),
    ],
)
def test_reader_rejects_duplicate_or_non_snp_loci(tmp_path, rows, match):
    path = tmp_path / "invalid.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=7,length=100>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\n"
        + "".join(
            f"7\t{pos}\t.\tA\t{alt}\t.\tPASS\t.\tGT\t0|1\t0|0\n" for pos, alt in rows
        )
    )
    with pytest.raises(ValueError, match=match):
        read_phased_biallelic_vcf(path)


def test_development_summary_recomputes_relative_coordinate_errors():
    records = [
        {
            "seed": i,
            "coordinate_origin": 10_000,
            "sequence_length": 1000,
            "source_record_sha256": "a" * 64,
            "source_tree_sha256": "b" * 64,
            "truth": {"position": 10_200},
            "predictions": {
                "a": {
                    "position": 250,
                    "window_snps": 50,
                    "minimum_minor_frequency": 0,
                    "informative": True,
                    "maximum_ties": 1,
                }
            },
            "absolute_errors": {"a": 50},
        }
        for i in range(48)
    ]
    assert (
        summarize_development(records)["development_metrics"]["a"][
            "median_absolute_error"
        ]
        == 50
    )
    records[0]["absolute_errors"]["a"] = -1
    with pytest.raises(ValueError, match="finite and valid"):
        summarize_development(records)
    records[0]["absolute_errors"]["a"] = 0
    with pytest.raises(ValueError, match="differ from predictions"):
        summarize_development(records)
    records[0]["absolute_errors"]["a"] = 50
    records[0]["predictions"]["a"]["window_snps"] = 100
    with pytest.raises(ValueError, match="configurations differ"):
        summarize_development(records)


def test_summary_requires_reference_raw_report_agreement(tmp_path):
    case = dict(
        case=0,
        seed=1,
        history="constant",
        selection_coefficient=0.01,
        origin_generations=1000,
        position=160_000,
    )
    design = tmp_path / "design.json"
    design.write_text(json.dumps({"cases": [case], "predictor": {"method": "test"}}))
    folder = tmp_path / "case_000"
    folder.mkdir()

    def digest(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()

    (folder / "status.json").write_text(
        json.dumps(
            {
                "case": 0,
                "seed": 1,
                "state": "COMPLETED",
                "design_sha256": digest(design),
            }
        )
    )
    (folder / "prediction.json").write_text(
        json.dumps(
            {
                "method": "test",
                "localization": {"truth_used_for_inference": False},
                "interval": [100_000, 220_000],
                "center": 160_000,
                "consensus_eligible": True,
                "contig_length": 1_000_000,
                "accepted_biallelic_sites": 100,
            }
        )
    )
    (folder / "candidate.bed").write_text("7\t100000\t220000\n")
    (folder / "RAiSD_Report.reference").write_text("900001 0 0 0 0 0 1\n")
    (folder / "RAiSD_Info.reference").write_text("| snps 100\n")
    record = {
        **case,
        "method": "test",
        "absolute_error": 0,
        "reference_position": 160_000,
        "reference_absolute_error": 0,
        "covered": True,
        "interval_width": 120_000,
        "consensus_eligible": True,
        "prediction_sha256": digest(folder / "prediction.json"),
        "bed_sha256": digest(folder / "candidate.bed"),
    }
    metrics_path = folder / "metrics.json"
    metrics_path.write_text(json.dumps(record))
    with pytest.raises(ValueError, match="Stored metrics differ"):
        summarize(design, tmp_path)
    record.update(reference_position=900_000, reference_absolute_error=740_000)
    metrics_path.write_text(json.dumps(record))
    result = summarize(design, tmp_path)
    assert (
        result["completed_only_metrics"]["reference_median_absolute_error"] == 740_000
    )
    record["reference_report_sha256"] = "0" * 64
    metrics_path.write_text(json.dumps(record))
    with pytest.raises(ValueError, match="Reference artifact hash"):
        summarize(design, tmp_path)
    with pytest.raises(ValueError, match="finite in-contig"):
        highest_peak_position(np.array([[1, 0, 0, 0, 0, 0, np.nan]]), 100)
