from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path

from ghist_relatedness.validation import validate_simulation_splits

import joblib
import numpy as np

try:
    from scripts.evaluate_human_background_robust import FINAL_REGIMES
    from scripts.fit_human_background_corrected import (
        correct_records,
        metrics,
        blended_probabilities,
    )
    from scripts.fit_human_background_deconvolution import current_baseline_metrics
    from scripts.fit_human_background_robust import (
        augment_records,
    )
    from scripts.validate_human_background_correction import read_pedsim_with_lengths
except ModuleNotFoundError:
    from evaluate_human_background_robust import FINAL_REGIMES
    from fit_human_background_corrected import (
        correct_records,
        metrics,
        blended_probabilities,
    )
    from fit_human_background_deconvolution import current_baseline_metrics
    from fit_human_background_robust import augment_records
    from validate_human_background_correction import read_pedsim_with_lengths


def aggregate(runs: list[dict[str, object]]) -> dict[str, float]:
    exact = np.asarray([float(row["exact_accuracy"]) for row in runs])
    degree = np.asarray([float(row["meioses_accuracy"]) for row in runs])
    return {
        "mean_exact_accuracy": float(exact.mean()),
        "minimum_exact_accuracy": float(exact.min()),
        "mean_meioses_accuracy": float(degree.mean()),
        "minimum_meioses_accuracy": float(degree.min()),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-model", type=Path, required=True)
    parser.add_argument("--baseline-model", type=Path, required=True)
    parser.add_argument("--holdouts", type=Path, nargs="+", required=True)
    parser.add_argument("--correction-reference", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--genome-cm", type=float, default=3346.298892538)
    parser.add_argument("--seed", type=int, default=121226)
    args = parser.parse_args()
    split_manifest = validate_simulation_splits(
        {"holdouts": args.holdouts, "correction_reference": [args.correction_reference]}
    )
    correction_reference, _ = read_pedsim_with_lengths(
        args.correction_reference, args.genome_cm
    )

    candidate = joblib.load(args.candidate_model)
    baseline = joblib.load(args.baseline_model)
    minimum_cm = min(
        int(name.rsplit("_", 1)[1])
        for name in candidate["feature_names"]
        if name.startswith("segments_ge_")
    )
    baseline_minimum_cm = min(
        int(name.rsplit("_", 1)[1])
        for name in baseline["feature_names"]
        if name.startswith("segments_ge_")
    )
    candidate_runs = []
    baseline_runs = []
    paired = []
    for replicate, path in enumerate(args.holdouts):
        reference, lengths = read_pedsim_with_lengths(path, args.genome_cm)
        for regime_index, regime in enumerate(FINAL_REGIMES):
            augmented = augment_records(
                reference,
                lengths,
                regime,
                np.random.default_rng(args.seed + 100 * replicate + regime_index),
                args.genome_cm,
            )
            corrected = correct_records(
                augmented, correction_reference, minimum_cm, args.genome_cm
            )
            classes, probabilities = blended_probabilities(
                candidate["corrected_payload"],
                baseline["model"],
                baseline["feature_names"],
                corrected,
                candidate["feature_names"],
                float(candidate["corrected_model_weight"]),
            )
            candidate_result = metrics(classes, probabilities, corrected)
            baseline_result = current_baseline_metrics(
                baseline["model"],
                baseline["feature_names"],
                augmented,
                correction_reference,
                baseline_minimum_cm,
                args.genome_cm,
            )
            identifiers = {
                "replicate": replicate,
                "source": str(path),
                "regime": asdict(regime),
            }
            candidate_runs.append({**identifiers, **candidate_result})
            baseline_runs.append({**identifiers, **baseline_result})
            paired.append(
                {
                    **identifiers,
                    "exact_difference": candidate_result["exact_accuracy"]
                    - baseline_result["exact_accuracy"],
                    "meioses_difference": candidate_result["meioses_accuracy"]
                    - baseline_result["meioses_accuracy"],
                }
            )

    exact_differences = np.asarray([float(row["exact_difference"]) for row in paired])
    degree_differences = np.asarray(
        [float(row["meioses_difference"]) for row in paired]
    )
    gate = {
        "mean_exact_gain_at_least_0.02": float(exact_differences.mean()) >= 0.02,
        "exact_gain_positive_in_at_least_80_percent": float(
            np.mean(exact_differences > 0)
        )
        >= 0.80,
        "worst_exact_regression_no_more_than_0.01": float(exact_differences.min())
        >= -0.01,
        "mean_meioses_regression_no_more_than_0.005": float(degree_differences.mean())
        >= -0.005,
    }
    gate["passed"] = all(gate.values())
    summary = {
        "state": "COMPLETED",
        "design": "preregistered paired comparison on newly generated Ped-sim replicates and unseen background regimes",
        "holdout_replicates": len(args.holdouts),
        "independent_correction_reference": True,
        "split_manifest": split_manifest,
        "regimes_share_pedigrees": True,
        "background_regimes": [asdict(regime) for regime in FINAL_REGIMES],
        "comparisons": len(paired),
        "candidate": aggregate(candidate_runs),
        "baseline": aggregate(baseline_runs),
        "mean_paired_exact_gain": float(exact_differences.mean()),
        "minimum_paired_exact_gain": float(exact_differences.min()),
        "positive_exact_comparisons": int(np.sum(exact_differences > 0)),
        "mean_paired_meioses_gain": float(degree_differences.mean()),
        "gate": gate,
        "candidate_runs": candidate_runs,
        "baseline_runs": baseline_runs,
        "paired_differences": paired,
        "submission_count": 0,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2) + "\n")
    print(
        json.dumps(
            {key: value for key, value in summary.items() if not key.endswith("runs")},
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
