"""Create an explicit versioned Task 9 candidate after fixed confirmation gates."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import math

import joblib
import sklearn

from ghist_relatedness import ibd
from ghist_relatedness.labels import IBD_LABEL_MEIOSES
from ghist_relatedness.task9_release import sha256, validate_payload


def release_gates(training, evidence, operational, model_hash):
    by_scenario = {result["scenario"]: result for result in evidence}
    if set(by_scenario) != {"confirmation_interference", "confirmation_poisson"} or len(evidence) != 2:
        raise ValueError("Require the two fixed untouched confirmation scenarios")
    matched = by_scenario["confirmation_interference"]["runs"]["unchanged"]
    train_hashes = {row["sha256"] for row in training["source_splits"]["calibration"]}
    evidence_hashes = {row["segment_sha256"] for row in evidence}
    primary = [row["runs"]["unchanged"] for row in evidence]
    expected_by_hash = {row["segment_sha256"]: row["runs"]["unchanged"] for row in evidence}
    return {
        # Preserve historical numerical gates; also require the shifted confirmation.
        "matched_confirmation_at_least_10000_pairs": matched["pairs"] >= 10000,
        "both_exact_accuracy_at_least_80pct": all(row["accuracy"] >= .80 for row in primary),
        "both_balanced_accuracy_at_least_75pct": all(row["balanced_accuracy"] >= .75 for row in primary),
        "both_meioses_accuracy_at_least_84pct": all(row["meioses_accuracy"] >= .84 for row in primary),
        "every_map_stress_has_probability_support": all(
            run["zero_accepted_probability_pairs"] == 0 for row in evidence for run in row["runs"].values()),
        "exactly_current_training_labels": set(training["labels"]) == set(IBD_LABEL_MEIOSES),
        "no_answer_or_confirmation_selection": training["public_answers_used"] is False
            and training["new_evaluation_used_for_selection"] is False
            and all(row["public_answers_used"] is False for row in evidence),
        "distinct_calibration_and_confirmation": len(evidence_hashes) == 2 and not train_hashes & evidence_hashes,
        "same_fixed_model": all(row["model_sha256"] == model_hash for row in evidence)
            and operational["model_sha256"] == model_hash,
        "production_pair_coverage_and_format": operational.get("passed") is True
            and len(operational["datasets"]) == 2
            and {row["source_sha256"] for row in operational["datasets"]} == evidence_hashes
            and all(row["permutation_invariant"] is True
                    and set(row["cases"]) == {"official_csv", "permuted_reversed_bom_csv"}
                    and all(case["validation"]["format_valid"] is True
                            and case["validation"]["input_pair_coverage_verified"] is True
                            and case["validation"]["pairs"] == row["expected_pairs"]
                            for case in row["cases"].values()) for row in operational["datasets"]),
        "production_matches_fixed_confirmation": len(operational["datasets"]) == 2 and all(
            row["source_sha256"] in expected_by_hash
            and row["expected_pairs"] == expected_by_hash[row["source_sha256"]]["pairs"]
            and all(math.isclose(case[metric], expected_by_hash[row["source_sha256"]][metric], abs_tol=1e-12)
                    for case in row["cases"].values() for metric in ("accuracy", "meioses_accuracy"))
            for row in operational["datasets"]),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--training", type=Path, required=True)
    parser.add_argument("--confirmation", type=Path, nargs=2, required=True)
    parser.add_argument("--operational-report", type=Path, required=True)
    parser.add_argument("--sex-map", type=Path, required=True)
    parser.add_argument("--artifact-root", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--release-id", default="2026-11-20-ten-labels")
    args = parser.parse_args()
    model_hash = sha256(args.model)
    training = json.loads(args.training.read_text())
    evidence = [json.loads(path.read_text()) for path in args.confirmation]
    operational = json.loads(args.operational_report.read_text())
    gates = release_gates(training, evidence, operational, model_hash)
    if not all(gates.values()):
        raise ValueError(f"Task 9 release gates failed: {gates}")
    payload = joblib.load(args.model)
    model, features = validate_payload(payload)
    if set(model.classes_) != set(IBD_LABEL_MEIOSES):
        raise ValueError("Release model has retired or missing labels")
    destination = args.artifact_root / "artifacts/task9_releases" / args.release_id / "model.joblib"
    if destination.exists() and sha256(destination) != model_hash:
        raise ValueError("Refusing to overwrite a different versioned model")
    destination.parent.mkdir(parents=True, exist_ok=True)
    if not destination.exists():
        shutil.copyfile(args.model, destination)
    manifest = {
        "schema": "ghist.task9.release.v1", "release_id": args.release_id,
        "status": "VALIDATED_RELEASE_CANDIDATE", "default_model_changed": False,
        "pair_domain": "unordered_pairs_present_in_segment_csv",
        "model": {"path": str(destination.relative_to(args.artifact_root)), "sha256": model_hash},
        "accepted_labels": sorted(IBD_LABEL_MEIOSES), "features": list(features),
        "genome_cm": 3350.933402, "sex_map_sha256": sha256(args.sex_map),
        "feature_contract_sha256": {name: sha256(Path(ibd.__file__).with_name(name))
                                    for name in ("ibd.py", "labels.py")},
        "environment": {"scikit_learn": sklearn.__version__}, "gates": gates,
        "evidence_sha256": {str(path): sha256(path) for path in
                            (args.training, *args.confirmation, args.operational_report)},
        "validation_scope": "Exact Ped-sim segment CSV; independent pedigrees; fixed matched and shifted crossover models",
        "limitations": ["Class probabilities are not calibrated posteriors",
                        "Distant-cousin labels remain difficult; accuracy is conditional on input pairs with IBD",
                        "No current challenge-data inference or submission performed"],
        "historical_freeze_preserved": True, "submission_count": 0,
    }
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    if args.manifest.exists() and json.loads(args.manifest.read_text()) != manifest:
        raise ValueError("Refusing to overwrite a different release manifest")
    args.manifest.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
