from __future__ import annotations

import argparse
import json
from pathlib import Path

import joblib
import numpy as np

try:
    from scripts.fit_task9_degree_factorized import aligned, combined_probabilities, metrics
    from scripts.fit_task9_ibd_classifier import load_simulations
except ModuleNotFoundError:
    from fit_task9_degree_factorized import aligned, combined_probabilities, metrics
    from fit_task9_ibd_classifier import load_simulations
from ghist_relatedness.ibd import matrix, read_sex_map


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-model", type=Path, required=True)
    parser.add_argument("--baseline-model", type=Path, required=True)
    parser.add_argument("--holdouts", type=Path, nargs="+", required=True)
    parser.add_argument("--sex-map", type=Path, required=True)
    parser.add_argument("--genome-cm", type=float, default=3350.933402)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    candidate = joblib.load(args.candidate_model)
    baseline = joblib.load(args.baseline_model)["model"]
    sex_map = read_sex_map(args.sex_map)
    runs = []
    for path in args.holdouts:
        records = load_simulations([path], args.genome_cm, sex_map)
        values, _ = matrix(records, candidate["features"])
        truth = np.asarray([record["label"] for record in records])
        candidate_probabilities = combined_probabilities(
            candidate["bundle"], baseline, values, candidate["winner"]
        )
        baseline_probabilities = aligned(baseline.classes_, baseline.predict_proba(values))
        candidate_metrics = metrics(candidate_probabilities, truth)
        baseline_metrics = metrics(baseline_probabilities, truth)
        runs.append(
            {
                "source": str(path),
                "candidate": candidate_metrics,
                "baseline": baseline_metrics,
                "exact_difference": candidate_metrics["accuracy"] - baseline_metrics["accuracy"],
                "balanced_difference": candidate_metrics["balanced_accuracy"] - baseline_metrics["balanced_accuracy"],
                "meioses_difference": candidate_metrics["meioses_accuracy"] - baseline_metrics["meioses_accuracy"],
            }
        )
    exact = np.asarray([run["exact_difference"] for run in runs])
    balanced = np.asarray([run["balanced_difference"] for run in runs])
    degree = np.asarray([run["meioses_difference"] for run in runs])
    gate = {
        "mean_exact_gain_at_least_0.002": float(exact.mean()) >= 0.002,
        "at_least_three_of_four_exact_gains_positive": int(np.sum(exact > 0)) >= 3,
        "mean_balanced_accuracy_non_regression": float(balanced.mean()) >= 0,
        "mean_meioses_non_regression": float(degree.mean()) >= 0,
    }
    gate["passed"] = all(gate.values())
    summary = {
        "state": "COMPLETED",
        "holdout_replicates": len(runs),
        "candidate_mean_exact_accuracy": float(np.mean([run["candidate"]["accuracy"] for run in runs])),
        "baseline_mean_exact_accuracy": float(np.mean([run["baseline"]["accuracy"] for run in runs])),
        "mean_paired_exact_gain": float(exact.mean()),
        "minimum_paired_exact_gain": float(exact.min()),
        "positive_exact_replicates": int(np.sum(exact > 0)),
        "mean_paired_balanced_gain": float(balanced.mean()),
        "candidate_mean_meioses_accuracy": float(np.mean([run["candidate"]["meioses_accuracy"] for run in runs])),
        "baseline_mean_meioses_accuracy": float(np.mean([run["baseline"]["meioses_accuracy"] for run in runs])),
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
