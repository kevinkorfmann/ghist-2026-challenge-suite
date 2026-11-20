from __future__ import annotations

import gzip
import json

import pytest

from scripts.validate_task10_complete_pipeline import (
    family_truth,
    normalize_vcf,
    save_result,
    score_predictions,
)


def test_family_truth_uses_shared_genealogy_not_sample_labels():
    parents = {"grandpa": (), "grandma": (), "mate1": (), "mate2": (),
               "parent1": ("grandpa", "grandma"), "parent2": ("grandpa", "grandma"),
               "a": ("parent1", "mate1"), "b": ("parent1", "mate1"),
               "c": ("parent2", "mate2")}
    rows = family_truth(parents, {"a": "sample_0", "b": "sample_1", "c": "sample_2"})
    assert [row["relation"] for row in rows] == ["siblings", "full_1st_cousin", "full_1st_cousin"]


def test_complete_three_branch_family_truth_has_all_153_pairs():
    from collections import Counter

    parents = {"root1": (), "root2": ()}
    identifiers = {}
    for branch in range(3):
        child, mate = f"child{branch}", f"mate{branch}"
        parents[child] = ("root1", "root2")
        parents[mate] = ()
        for subbranch in range(3):
            grandchild, spouse = f"g{branch}_{subbranch}", f"s{branch}_{subbranch}"
            parents[grandchild] = (child, mate)
            parents[spouse] = ()
            for sibling in range(2):
                sample = f"leaf{branch}_{subbranch}_{sibling}"
                parents[sample] = (grandchild, spouse)
                identifiers[sample] = f"sample_{len(identifiers)}"
    truth = family_truth(parents, identifiers)
    assert len(truth) == 153
    assert Counter(row["relation"] for row in truth) == {
        "siblings": 9, "full_1st_cousin": 36, "full_2nd_cousin": 108,
    }


def test_missing_zero_ibd_pair_counts_as_wrong():
    truth = [{"indID1": "a", "indID2": "b", "relation": "siblings"},
             {"indID1": "a", "indID2": "c", "relation": "full_2nd_cousin"}]
    score = score_predictions(truth, truth[:1])
    assert score["missing_pairs"] == 1
    assert score["exact_accuracy"] == 0.5
    assert not score["all_pairs_present"]


def test_reversed_duplicate_and_unknown_pairs_are_rejected():
    row = {"indID1": "a", "indID2": "b", "relation": "siblings"}
    with pytest.raises(ValueError, match="duplicate"):
        score_predictions([row], [row, {**row, "indID1": "b", "indID2": "a"}])
    with pytest.raises(ValueError, match="unknown"):
        score_predictions([row], [{**row, "indID2": "c"}])


@pytest.mark.parametrize("call", ["0/1", "./.", "0|."])
def test_oracle_phase_contract_rejects_unphased_and_missing(tmp_path, call):
    source = tmp_path / "source.vcf.gz"
    with gzip.open(source, "wt") as handle:
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\n")
        handle.write(f"1\t1\t.\tA\tC\t.\tPASS\t.\tGT\t{call}\n")
    with pytest.raises(ValueError, match="rejects genotype"):
        normalize_vcf(source, tmp_path / "out.vcf.gz")


def test_result_binds_exact_generated_truth_inputs_and_candidate(tmp_path):
    from ghist_relatedness.provenance import sha256, verify_hashes

    (tmp_path / "preregistered_design.json").write_text('{"seed": 17}\n')
    (tmp_path / "truth.json").write_text('[]\n')
    (tmp_path / "prepared.vcf.gz").write_bytes(b"synthetic genotype bytes")
    (tmp_path / "candidate.csv").write_text("indID1,indID2,relation\na,b,siblings\n")
    save_result(tmp_path, {"state": "COMPLETED"})
    record = json.loads((tmp_path / "result.json").read_text())
    assert record["protocol_sha256"] == sha256(tmp_path / "preregistered_design.json")
    assert set(record["output_hashes"]) == {
        "preregistered_design.json", "truth.json", "prepared.vcf.gz", "candidate.csv"
    }
    bindings = {str(tmp_path / name): digest for name, digest in record["output_hashes"].items()}
    verify_hashes(bindings, context="Retained Task 10 evidence")
    (tmp_path / "candidate.csv").write_text("indID1,indID2,relation\na,b,half_sibling\n")
    with pytest.raises(ValueError, match="SHA-256 mismatch.*candidate.csv"):
        verify_hashes(bindings, context="Retained Task 10 evidence")
