from __future__ import annotations

import argparse
import copy
import json
from dataclasses import asdict
from pathlib import Path

from ghist_relatedness.validation import validate_simulation_splits

import joblib
import numpy as np

try:
    from scripts.fit_human_background_robust import (
        BackgroundRegime,
        augment_records,
        metrics,
    )
    from scripts.fit_human_family_relationships import background_correct
    from scripts.validate_human_background_correction import (
        accuracy,
        read_pedsim_with_lengths,
        trim_shortest_background_shape,
    )
except ModuleNotFoundError:
    from fit_human_background_robust import BackgroundRegime, augment_records, metrics
    from fit_human_family_relationships import background_correct
    from validate_human_background_correction import (
        accuracy,
        read_pedsim_with_lengths,
        trim_shortest_background_shape,
    )


FINAL_REGIMES = (
    BackgroundRegime(1.50, 20.0, 0.28, 0.10, 3.0),
    BackgroundRegime(2.75, 8.0, 0.40, 0.18, 2.0),
    BackgroundRegime(4.25, 25.0, 0.18, 0.04, 3.5),
    BackgroundRegime(5.50, 12.0, 0.52, 0.15, 2.8),
    BackgroundRegime(7.00, 18.0, 0.35, 0.22, 2.2),
)


def summarize_runs(runs: list[dict[str, object]]) -> dict[str, float]:
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
    parser.add_argument("--robust-model", type=Path, required=True)
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

    robust = joblib.load(args.robust_model)
    baseline = joblib.load(args.baseline_model)
    robust_runs = []
    baseline_runs = []
    paired_differences = []
    for replicate, path in enumerate(args.holdouts):
        reference, source_lengths = read_pedsim_with_lengths(path, args.genome_cm)
        minimum_cm = min(
            int(name.rsplit("_", 1)[1])
            for name in baseline["feature_names"]
            if name.startswith("segments_ge_")
        )
        for regime_index, regime in enumerate(FINAL_REGIMES):
            augmented = augment_records(
                reference,
                source_lengths,
                regime,
                np.random.default_rng(args.seed + 100 * replicate + regime_index),
                args.genome_cm,
            )
            robust_result = metrics(
                robust["payload"], augmented, robust["feature_names"]
            )
            corrected = copy.deepcopy(augmented)
            offsets = background_correct(
                corrected, correction_reference, quantile=0.25, genome_cm=args.genome_cm
            )
            trim_shortest_background_shape(
                corrected,
                float(offsets[f"segments_ge_{minimum_cm}"]),
                minimum_cm,
            )
            baseline_exact, baseline_degree = accuracy(
                baseline["model"], corrected, baseline["feature_names"]
            )
            identifiers = {
                "replicate": replicate,
                "source": str(path),
                "regime": asdict(regime),
            }
            robust_runs.append({**identifiers, **robust_result})
            baseline_runs.append(
                {
                    **identifiers,
                    "exact_accuracy": baseline_exact,
                    "meioses_accuracy": baseline_degree,
                }
            )
            paired_differences.append(
                {
                    **identifiers,
                    "exact_difference": robust_result["exact_accuracy"]
                    - baseline_exact,
                    "meioses_difference": robust_result["meioses_accuracy"]
                    - baseline_degree,
                }
            )

    exact_differences = np.asarray(
        [float(row["exact_difference"]) for row in paired_differences]
    )
    degree_differences = np.asarray(
        [float(row["meioses_difference"]) for row in paired_differences]
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
        "comparisons": len(paired_differences),
        "robust": summarize_runs(robust_runs),
        "baseline": summarize_runs(baseline_runs),
        "mean_paired_exact_gain": float(exact_differences.mean()),
        "minimum_paired_exact_gain": float(exact_differences.min()),
        "mean_paired_meioses_gain": float(degree_differences.mean()),
        "gate": gate,
        "robust_runs": robust_runs,
        "baseline_runs": baseline_runs,
        "paired_differences": paired_differences,
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
