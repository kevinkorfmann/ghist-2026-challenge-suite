"""Reproduce retrospective development replay; never refit or read challenge data."""

from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import shutil

from scripts.consensus_task11_dog_pedigrees import (
    MEIOSES,
    bonsai_label,
    build_consensus,
    pair_key,
)
from scripts.consensus_task11_guarded_candidate import (
    build_guarded_candidate,
    source_hashes,
)


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
HISTORICAL = ROOT / "experiments/complete_pipeline_2026_11_20/task11"
SOURCES = (
    "scripts/consensus_task11_dog_pedigrees.py",
    "scripts/consensus_task11_guarded_candidate.py",
    "scripts/freeze_task11_multifamily_inbreeding.py",
    "tests/test_task11_consensus.py",
    "tests/test_task11_freeze_guard.py",
    "src/ghist_relatedness/labels.py",
)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")


def main():
    protocol = json.loads((HERE / "policy_protocol.json").read_text())
    original_protocol = json.loads(
        (HISTORICAL / "results/seed_2026112711/protocol.json").read_text()
    )
    original_sources = {
        str(path.relative_to(HERE / "original_source")): sha(path)
        for path in (HERE / "original_source").rglob("*.py")
    }
    assert all(
        original_protocol["sources"][relative] == value
        for relative, value in original_sources.items()
    )
    retained_hashes = {
        str(path.relative_to(ROOT)): sha(path)
        for path in sorted(HISTORICAL.rglob("*"))
        if path.is_file()
    }
    verified_retained_outputs = 0
    cohorts = []
    differences = []
    for seed in protocol["development_cohort_seeds"]:
        source = HISTORICAL / "results" / f"seed_{seed}"
        result = json.loads((source / "result.json").read_text())
        for relative, expected in result["output_hashes"].items():
            path = source / relative
            if path.is_file():
                assert sha(path) == expected, path
                verified_retained_outputs += 1
        original = json.loads((source / "consensus/predictions.json").read_text())
        original_by_pair = {pair_key(row): row for row in original}
        forest = [
            json.loads((source / f"pair_predictions_hbd{threshold}.json").read_text())
            for threshold in (2, 5, 10)
        ]
        topology = json.loads((source / "bonsai.json").read_text())
        forest_by_pair = [{pair_key(row): row for row in run} for run in forest]
        label_support_audit = []
        for old in original:
            pair = pair_key(old)
            label_support_audit.append(
                {
                    "indID1": pair[0],
                    "indID2": pair[1],
                    "relation": old["relation"],
                    "source": old["source"],
                    "historical_confidence": old["confidence"],
                    "confidence": sum(
                        run[pair]["probabilities"][old["relation"]]
                        for run in forest_by_pair
                    )
                    / len(forest),
                    "scope": "Metadata-only audit of unchanged historical label; not a new prediction.",
                }
            )
        rejection = None
        try:
            audited, audited_summary = build_consensus(forest, topology)
        except ValueError as error:
            if str(error) != "No probability support for topology-constrained labels":
                raise
            audited, audited_summary = None, None
            rejection = str(error)
        zero_constraints = []
        for row in topology["relationships"]:
            label, meioses = bonsai_label(row)
            if label is not None or meioses is None:
                continue
            pair = pair_key(row)
            for threshold, run in zip((2, 5, 10), forest_by_pair, strict=True):
                if (
                    sum(
                        value
                        for name, value in run[pair]["probabilities"].items()
                        if MEIOSES[name] == meioses
                    )
                    == 0
                ):
                    zero_constraints.append(
                        {
                            "pair": pair,
                            "hbd_threshold_cm": threshold,
                            "topology_meioses": meioses,
                        }
                    )
        guarded, guarded_summary = build_guarded_candidate(forest, topology)
        guarded_summary["source_sha256"] = source_hashes()
        guarded_summary["input_sha256"] = {
            str(path.relative_to(ROOT)): sha(path)
            for path in [
                source / f"pair_predictions_hbd{threshold}.json"
                for threshold in (2, 5, 10)
            ]
            + [source / "bonsai.json"]
        }
        destination = HERE / f"seed_{seed}"
        if audited is not None:
            assert all(
                row["relation"] == original_by_pair[pair_key(row)]["relation"]
                for row in audited
            )
            assert all(
                row["confidence"] == row["probabilities"][row["relation"]]
                for row in audited
            )
            write(destination / "legacy_audited/predictions.json", audited)
            write(destination / "legacy_audited/summary.json", audited_summary)
        else:
            write(
                destination / "legacy_audited/rejected.json",
                {
                    "reason": rejection,
                    "zero_probability_constraints": zero_constraints,
                    "candidate_exported": False,
                },
            )
        write(destination / "original_label_support_audit.json", label_support_audit)
        write(destination / "guarded_development/development_predictions.json", guarded)
        write(
            destination / "guarded_development/development_summary.json",
            guarded_summary,
        )
        for row in label_support_audit:
            old = original_by_pair[pair_key(row)]
            if abs(row["confidence"] - old["confidence"]) > 1e-12:
                differences.append(
                    {
                        "seed": seed,
                        "indID1": row["indID1"],
                        "indID2": row["indID2"],
                        "unchanged_relation": row["relation"],
                        "source": row["source"],
                        "old_confidence": old["confidence"],
                        "corrected_final_label_probability": row["confidence"],
                    }
                )
        cohorts.append(
            {
                "seed": seed,
                "pairs": len(original),
                "original_consensus_correct": sum(
                    row["relation"] == "siblings" for row in original
                ),
                "audited_legacy_consensus_correct": sum(
                    row["relation"] == "siblings" for row in audited
                )
                if audited is not None
                else None,
                "legacy_guard_rejection": rejection,
                "zero_probability_constraints": zero_constraints,
                "legacy_label_changes_on_successful_replay": 0
                if audited is not None
                else None,
                "raw_forest_proposals_correct": sum(
                    row["proposed_relation"] == "siblings" for row in guarded
                ),
                "guarded_provisional_pairs": guarded_summary["provisional_pairs"],
                "guarded_abstained_pairs": guarded_summary["abstained_pairs"],
                "guarded_correct_counting_abstentions_as_wrong": sum(
                    row["relation"] == "siblings" for row in guarded
                ),
                "guarded_abstention_reasons": guarded_summary["abstention_reasons"],
                "raw_hbd_stability": sum(
                    len(set(row["raw_sensitivity_labels"])) == 1 for row in guarded
                )
                / len(guarded),
                "historical_topology_conditioned_stability": json.loads(
                    (source / "consensus/summary.json").read_text()
                )["hbd_sensitivity_stability"],
            }
        )
    assert all(sha(ROOT / path) == value for path, value in retained_hashes.items())
    with (HERE / "confidence_corrections.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(differences[0]))
        writer.writeheader()
        writer.writerows(differences)
    sources = {relative: sha(ROOT / relative) for relative in SOURCES}
    for relative in SOURCES:
        destination = HERE / "source" / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / relative, destination)
    summary = {
        "scope": "Retrospective development diagnostics on previously inspected sibling-only cohorts; no independent confirmation or accuracy certification.",
        "independent_families": 16,
        "development_cohorts": cohorts,
        "pairs": sum(row["pairs"] for row in cohorts),
        "legacy_label_changes_on_successful_replay": 0,
        "legacy_guard_rejected_cohorts": sum(
            row["legacy_guard_rejection"] is not None for row in cohorts
        ),
        "confidence_corrections": len(differences),
        "confidence_corrections_by_source": dict(
            Counter(row["source"] for row in differences)
        ),
        "original_smoke_output_hashes_verified": verified_retained_outputs,
        "historical_files_hashes_unchanged": len(retained_hashes),
        "historical_file_sha256": retained_hashes,
        "source_sha256": sources,
        "original_source_sha256": original_sources,
        "development_policy_protocol_sha256": sha(HERE / "policy_protocol.json"),
        "submission_ready": False,
        "submission_count": 0,
        "release_validation_evidence_created": False,
    }
    write(HERE / "audit_summary.json", summary)
    print(
        json.dumps(
            {
                key: value
                for key, value in summary.items()
                if not key.endswith("sha256")
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
