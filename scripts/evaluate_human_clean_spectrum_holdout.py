from __future__ import annotations

import argparse
import json
from pathlib import Path

import joblib
import numpy as np

try:
    from scripts.fit_human_background_corrected import align_probabilities
    from scripts.fit_human_clean_spectrum import metrics
    from scripts.fit_human_family_relationships import LABELS, matrix, read_pedsim
except ModuleNotFoundError:
    from fit_human_background_corrected import align_probabilities
    from fit_human_clean_spectrum import metrics
    from fit_human_family_relationships import LABELS, matrix, read_pedsim


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-model", type=Path, required=True)
    parser.add_argument("--baseline-model", type=Path, required=True)
    parser.add_argument("--holdouts", type=Path, nargs="+", required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--genome-cm", type=float, default=3346.298892538)
    args = parser.parse_args()

    candidate = joblib.load(args.candidate_model)
    baseline = joblib.load(args.baseline_model)
    classes = np.asarray(LABELS)
    runs = []
    for replicate, path in enumerate(args.holdouts):
        records = read_pedsim(path, args.genome_cm)
        candidate_probabilities = candidate["candidate_model"].predict_proba(
            matrix(records, candidate["candidate_feature_names"])
        )
        baseline_probabilities = baseline["model"].predict_proba(
            matrix(records, baseline["feature_names"])
        )
        candidate_aligned = align_probabilities(
            np.asarray(candidate["candidate_model"].classes_),
            candidate_probabilities,
            classes,
        )
        baseline_aligned = align_probabilities(
            np.asarray(baseline["model"].classes_), baseline_probabilities, classes
        )
        weight = float(candidate["candidate_weight"])
        blended = weight * candidate_aligned + (1.0 - weight) * baseline_aligned
        candidate_result = metrics(classes, blended, records)
        baseline_result = metrics(classes, baseline_aligned, records)
        runs.append(
            {
                "replicate": replicate,
                "source": str(path),
                "candidate": candidate_result,
                "baseline": baseline_result,
                "exact_difference": candidate_result["exact_accuracy"]
                - baseline_result["exact_accuracy"],
                "meioses_difference": candidate_result["meioses_accuracy"]
                - baseline_result["meioses_accuracy"],
            }
        )
    exact = np.asarray([row["exact_difference"] for row in runs])
    degree = np.asarray([row["meioses_difference"] for row in runs])
    gate = {
        "mean_exact_gain_at_least_0.01": float(exact.mean()) >= 0.01,
        "exact_gain_positive_in_at_least_three_of_four": int(np.sum(exact > 0)) >= 3,
        "worst_exact_regression_no_more_than_0.005": float(exact.min()) >= -0.005,
        "mean_meioses_regression_no_more_than_0.005": float(degree.mean()) >= -0.005,
    }
    gate["passed"] = all(gate.values())
    summary = {
        "state": "COMPLETED",
        "design": "paired comparison on four newly generated untouched Ped-sim replicates",
        "holdout_replicates": len(runs),
        "candidate_mean_exact_accuracy": float(
            np.mean([row["candidate"]["exact_accuracy"] for row in runs])
        ),
        "baseline_mean_exact_accuracy": float(
            np.mean([row["baseline"]["exact_accuracy"] for row in runs])
        ),
        "mean_paired_exact_gain": float(exact.mean()),
        "minimum_paired_exact_gain": float(exact.min()),
        "positive_exact_replicates": int(np.sum(exact > 0)),
        "candidate_mean_meioses_accuracy": float(
            np.mean([row["candidate"]["meioses_accuracy"] for row in runs])
        ),
        "baseline_mean_meioses_accuracy": float(
            np.mean([row["baseline"]["meioses_accuracy"] for row in runs])
        ),
        "mean_paired_meioses_gain": float(degree.mean()),
        "gate": gate,
        "runs": runs,
        "submission_count": 0,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({key: value for key, value in summary.items() if key != "runs"}, indent=2))


if __name__ == "__main__":
    main()
