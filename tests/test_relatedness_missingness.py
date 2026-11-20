import gzip
import argparse
import json
from itertools import combinations

import numpy as np
import pytest

from scripts.prepare_hapibd_vcf import prepare_vcf
from scripts.validate_relatedness_missingness import mask_genotypes
from scripts import validate_relatedness_missingness as harness
from scripts.validate_task11_complete_pipeline import generate_cohort


def fixture(path):
    text = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\tc\n"
        + "".join(
            f"1\t{i}\t.\tA\tC\t.\tPASS\t.\tGT\t0|1\t1|0\t1|1\n" for i in range(1, 101)
        )
    )
    with gzip.open(path, "wt") as handle:
        handle.write(text)
    return text


def read(path):
    with gzip.open(path, "rt") as handle:
        return handle.read()


def test_masking_changes_only_selected_genotypes_and_is_reproducible(tmp_path):
    source, first, second = [tmp_path / name for name in ("src.gz", "a.gz", "b.gz")]
    original = fixture(source)
    a = mask_genotypes(source, first, 0.2, 20261120)
    b = mask_genotypes(source, second, 0.2, 20261120)
    assert a == b
    assert read(first) == read(second)
    assert 0 < a["masked_genotypes"] < 300
    observed_missing = 0
    for before, after in zip(
        original.splitlines(), read(first).splitlines(), strict=True
    ):
        if before.startswith("#"):
            assert before == after
            continue
        fields1, fields2 = before.split("\t"), after.split("\t")
        assert fields1[:9] == fields2[:9]
        for original_call, masked_call in zip(fields1[9:], fields2[9:], strict=True):
            assert masked_call in {original_call, ".|."}
            observed_missing += masked_call == ".|."
    assert observed_missing == a["masked_genotypes"]


def test_masking_zero_is_exact_control_and_policies_share_same_observations(tmp_path):
    source, zero, masked = [
        tmp_path / name for name in ("src.gz", "zero.gz", "masked.gz")
    ]
    original = fixture(source)
    assert mask_genotypes(source, zero, 0, 1)["masked_genotypes"] == 0
    assert read(zero) == original
    counts = mask_genotypes(source, masked, 0.2, 2)
    drop = prepare_vcf(masked, tmp_path / "drop.gz", missing_policy="drop-site")
    major = prepare_vcf(masked, tmp_path / "major.gz", missing_policy="major")
    assert drop["missing_records_dropped"] == counts["sites_with_missing_genotypes"]
    assert drop["imputed_alleles"] == 0
    assert major["records_written"] >= drop["records_written"]
    assert major["imputed_alleles"] > 0


def test_existing_missing_or_unphased_source_is_not_an_eligible_control(tmp_path):
    source, target = tmp_path / "src.gz", tmp_path / "out.gz"
    original = fixture(source)
    with gzip.open(source, "wt") as handle:
        handle.write(original.replace("0|1", "0/1", 1))
    with pytest.raises(ValueError, match="complete, phased"):
        mask_genotypes(source, target, 0.2, 1)


def baseline_fixture(tmp_path):
    baseline = tmp_path / "baseline"
    baseline.mkdir()
    map_path = tmp_path / "map"
    map_path.write_text("1 . 0 1\n1 . 0.2 10001\n")
    truth = generate_cohort(map_path, baseline, 2026112711)
    artifacts = {"map": map_path}
    for name in ("hapibd", "model"):
        artifacts[name] = tmp_path / name
        artifacts[name].write_bytes(name.encode())
    protocol = {
        "seed": truth["seed"],
        "fixed_seeds": [2026112711, 2026112811],
        "provenance": {
            name: {"path": str(path), "sha256": harness.sha256(path)}
            for name, path in artifacts.items()
        },
    }
    harness.write_json(baseline / "protocol.json", protocol)
    harness.write_json(baseline / "result.json", {
        "protocol_sha256": harness.sha256(baseline / "protocol.json"),
        "output_hashes": {
            name: harness.sha256(baseline / name)
            for name in ("protocol.json", "truth.json", "synthetic.vcf.gz")
        },
        "pair_model_by_hbd_threshold": {"5.0": {"correct": 48}},
        "family_partition": {"exactly_recovered_families": 8},
    })
    return baseline, truth, artifacts


@pytest.mark.parametrize("change", ["truth", "protocol", "result_protocol_digest", "missing_truth_digest"])
def test_baseline_identity_rejects_changed_truth_or_protocol_before_running(tmp_path, change):
    baseline, truth, _ = baseline_fixture(tmp_path)
    if change == "truth":
        truth["relationship"] = "half_sibling"
        harness.write_json(baseline / "truth.json", truth)
    elif change == "protocol":
        path = baseline / "protocol.json"
        record = json.loads(path.read_text())
        record["changed_after_baseline"] = True
        harness.write_json(path, record)
    else:
        path = baseline / "result.json"
        record = json.loads(path.read_text())
        if change == "result_protocol_digest":
            record["protocol_sha256"] = "0" * 64
        else:
            del record["output_hashes"]["truth.json"]
        harness.write_json(path, record)
    out = tmp_path / "out"
    with pytest.raises(ValueError, match="Baseline .*provenance|Baseline result/protocol linkage"):
        harness.run(argparse.Namespace(baseline=baseline, out=out, rate=0.005))
    assert not out.exists()


@pytest.mark.parametrize("change_model", [False, True])
def test_missingness_records_exact_outputs_and_rejects_inputs_changed_during_run(tmp_path, monkeypatch, change_model):
    baseline, truth, artifacts = baseline_fixture(tmp_path)
    pairs = [pair for family in truth["families"] for pair in combinations(family, 2)]

    class Frozen:
        def predict(self, features):
            return np.repeat("siblings", len(features))

    def detected_stage(name, command, out, stages):
        stages.append({"stage": name, "exit_code": 0})
        if change_model:
            artifacts["model"].write_bytes(b"replaced model")

    def cluster_fixture(pair_path, samples, out):
        clusters = {"families": truth["families"], "minimum_adjusted_rand": 1.0}
        harness.write_json(out / "clusters.json", clusters)
        return clusters

    monkeypatch.setattr(harness.joblib, "load", lambda path: Frozen())
    monkeypatch.setattr(harness, "run_command", detected_stage)
    monkeypatch.setattr(harness, "cluster_cohort", cluster_fixture)
    monkeypatch.setattr(harness, "read_observed", lambda *args: (np.zeros((48, 1)), pairs, None))
    out = tmp_path / "out"
    result = harness.run(argparse.Namespace(baseline=baseline, out=out, rate=0.005))
    for policy, record in result["policies"].items():
        assert record["status"] == ("FAILED" if change_model else "COMPLETED")
        if change_model:
            assert "SHA-256 mismatch" in record["failure"]["message"]
        else:
            assert record["input_provenance_verified"]
            assert record["pair_model"]["truth_pairs"] == 48
        for filename in ("prepared.vcf.gz", "pair_predictions.json", "clusters.json"):
            relative = f"{policy}/{filename}"
            assert result["output_hashes"][relative] == harness.sha256(out / relative)
    protocol = json.loads((out / "protocol.json").read_text())
    assert protocol["baseline_truth_sha256"] == harness.sha256(baseline / "truth.json")
    assert result["output_hashes"]["masked.vcf.gz"] == harness.sha256(out / "masked.vcf.gz")
