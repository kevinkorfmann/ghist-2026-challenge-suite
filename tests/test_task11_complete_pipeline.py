import argparse
import json
from itertools import combinations

import numpy as np
import pytest

from scripts import validate_task11_complete_pipeline as harness
from scripts.prepare_hapibd_vcf import prepare_vcf
from scripts.validate_task11_complete_pipeline import (
    family_ancestry,
    gamete,
    partition_metrics,
    relationship_metrics,
)


def test_joint_family_founder_haplotypes_are_disjoint_and_inheritance_is_valid():
    rng = np.random.default_rng(19011)
    cm = np.linspace(0, 100, 1001)
    previous_ids = set()
    offset = 0
    for mating in ("unrelated", "first_cousins", "siblings"):
        family, next_offset = family_ancestry(mating, cm, rng, offset)
        assert family.shape == (4, 2, 1001)
        assert family.min() >= offset
        assert family.max() < next_offset
        assert not (set(family.ravel()) & previous_ids)
        previous_ids.update(family.ravel())
        if mating == "unrelated":
            assert not np.any(family[:, 0] == family[:, 1])
        offset = next_offset


def test_gamete_preserves_parental_ancestry_at_each_marker():
    parent = np.array([[10, 11, 12], [20, 21, 22]])
    for seed in range(12):
        result = gamete(
            parent, np.array([0.0, 50.0, 100.0]), np.random.default_rng(seed)
        )
        assert np.all((result == parent[0]) | (result == parent[1]))


def test_partition_counts_missing_and_merged_pairs_instead_of_dropping_them():
    result = partition_metrics([["a", "b"], ["c", "d"]], [["a", "c"], ["b"]])
    assert result["missing_samples"] == ["d"]
    assert result["missed_within_pairs"] == 2
    assert result["false_cross_family_pairs"] == 1
    assert result["exactly_recovered_families"] == 0
    assert result["split_families"] == 1
    assert result["merged_predicted_families"] == 1


def test_partition_rejects_duplicate_or_unknown_samples():
    with pytest.raises(ValueError, match="Duplicate"):
        partition_metrics([["a", "b"]], [["a"], ["a", "b"]])
    with pytest.raises(ValueError, match="Unknown"):
        partition_metrics([["a", "b"]], [["a", "x"]])


def test_relationship_denominator_includes_missing_true_pairs_and_cross_family_errors():
    truth = {
        "families": [["a", "b"], ["c", "d"]],
        "samples": ["a", "b", "c", "d"],
        "family_mating_types": ["unrelated", "siblings"],
        "relationship": "siblings",
    }
    result = relationship_metrics(
        truth,
        [
            {"indID1": "b", "indID2": "a", "relation": "siblings"},
            {"indID1": "a", "indID2": "c", "relation": "siblings"},
        ],
    )
    assert result["truth_pairs"] == 2
    assert result["accuracy_counting_missing_as_wrong"] == 0.5
    assert result["missing_true_pairs"] == 1
    assert result["extra_cross_family_predictions"] == 1
    assert (
        result["by_parent_mating_type"]["siblings"][
            "accuracy_counting_missing_as_wrong"
        ]
        == 0
    )


def test_relationship_rejects_duplicate_predictions():
    truth = {
        "families": [["a", "b"]],
        "samples": ["a", "b"],
        "family_mating_types": ["unrelated"],
        "relationship": "siblings",
    }
    row = {"indID1": "a", "indID2": "b", "relation": "siblings"}
    with pytest.raises(ValueError, match="Duplicate prediction"):
        relationship_metrics(truth, [row, row])


@pytest.mark.parametrize("changed_input", [None, "model", "truth.json"])
def test_run_checks_original_inputs_before_accepting_completed_stages(tmp_path, monkeypatch, changed_input):
    map_path = tmp_path / "map"
    map_path.write_text("1 . 0 1\n1 . 0.2 10001\n")
    model_path, jar = tmp_path / "model", tmp_path / "hapibd.jar"
    model_path.write_bytes(b"frozen model")
    jar.write_bytes(b"detector")
    out = tmp_path / "out"

    def current_truth():
        return json.loads((out / "truth.json").read_text())

    def expected_pairs():
        return [pair for family in current_truth()["families"] for pair in combinations(family, 2)]

    def stage_fixture(name, command, directory, stages):
        stages.append({"stage": name, "exit_code": 0})
        if name == "prepare":
            prepare_vcf(out / "synthetic.vcf.gz", out / "prepared.vcf.gz")
        elif name == "segment_summary":
            (out / "summary").mkdir()
            (out / "summary/hbd_features.csv").write_text("individual,inbreeding_fraction\n")
        elif name == "consensus":
            harness.write_json(out / "consensus/predictions.json", [
                {"indID1": a, "indID2": b, "relation": "siblings"}
                for a, b in expected_pairs()
            ])
            if changed_input == "model":
                model_path.write_bytes(b"different model after inference")
            elif changed_input == "truth.json":
                truth = current_truth()
                truth["relationship"] = "half_sibling"
                harness.write_json(out / "truth.json", truth)

    def cluster_fixture(pair_path, samples, directory):
        clusters = {
            "families": current_truth()["families"],
            "state": "COMPLETED", "minimum_adjusted_rand": 1.0,
        }
        harness.write_json(out / "clusters.json", clusters)
        return clusters

    class Frozen:
        classes_ = np.array(["siblings"])

        def predict_proba(self, features):
            return np.ones((len(features), 1))

    monkeypatch.setattr(harness, "run_command", stage_fixture)
    monkeypatch.setattr(harness, "cluster_cohort", cluster_fixture)
    monkeypatch.setattr(harness.joblib, "load", lambda path: Frozen())
    monkeypatch.setattr(harness, "read_observed", lambda *args: (
        np.zeros((48, len(harness.feature_names()))), expected_pairs(), [0] * 48
    ))
    result = harness.run_pipeline(argparse.Namespace(
        map=map_path, hapibd=jar, model=model_path, seed=harness.SEEDS[0], out=out
    ))
    assert result["operational_stages_completed"] == (changed_input is None)
    if changed_input is None:
        assert result["input_provenance_verified"]
        assert result["consensus"]["truth_pairs"] == 48
    else:
        assert "SHA-256 mismatch" in result["failure"]["message"]
        assert not result.get("input_provenance_verified", False)
    assert result["provenance_schema"] == "relatedness-pipeline-v2"
    assert result["output_hashes"]["consensus/predictions.json"] == harness.sha256(out / "consensus/predictions.json")
