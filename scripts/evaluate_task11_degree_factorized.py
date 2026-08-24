from __future__ import annotations

import argparse
import json
from pathlib import Path

import joblib
import numpy as np

try:
    from scripts.fit_task11_degree_factorized import (
        aligned,
        degree_guarded_probabilities,
        factorized_probabilities,
        metrics,
        read_pedsim_rich,
    )
    from scripts.fit_task11_dog_relationships import feature_names
except ModuleNotFoundError:
    from fit_task11_degree_factorized import (
        aligned,
        degree_guarded_probabilities,
        factorized_probabilities,
        metrics,
        read_pedsim_rich,
    )
    from fit_task11_dog_relationships import feature_names


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-model", type=Path, required=True)
    parser.add_argument("--baseline-model", type=Path, required=True)
    parser.add_argument("--holdouts", type=Path, nargs="+", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    candidate = joblib.load(args.candidate_model)
    baseline = joblib.load(args.baseline_model)
    winner = candidate["winner"]
    bundle = candidate["bundle"]
    runs = []
    for path in args.holdouts:
        values, truth = read_pedsim_rich(path)
        factorized = factorized_probabilities(bundle, values)
        direct = aligned(
            np.asarray(bundle["exact"].classes_),
            bundle["exact"].predict_proba(values),
        )
        baseline_probabilities = aligned(
            np.asarray(baseline.classes_),
            baseline.predict_proba(values[:, : len(feature_names())]),
        )
        structured = winner["factor_weight"] * factorized + (1.0 - winner["factor_weight"]) * direct
        probabilities = winner["candidate_weight"] * structured + (1.0 - winner["candidate_weight"]) * baseline_probabilities
        probabilities = degree_guarded_probabilities(probabilities, baseline_probabilities)
        candidate_metrics = metrics(probabilities, truth)
        baseline_metrics = metrics(baseline_probabilities, truth)
        runs.append(
            {
                "source": str(path),
                "candidate": candidate_metrics,
                "baseline": baseline_metrics,
                "exact_difference": candidate_metrics["accuracy"] - baseline_metrics["accuracy"],
                "meioses_difference": candidate_metrics["meioses_accuracy"] - baseline_metrics["meioses_accuracy"],
            }
        )
    exact = np.asarray([run["exact_difference"] for run in runs])
    degree = np.asarray([run["meioses_difference"] for run in runs])
    gate = {
        "mean_exact_gain_at_least_0.005": float(exact.mean()) >= 0.005,
        "all_exact_gains_positive": bool(np.all(exact > 0)),
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
