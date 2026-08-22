from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path

import numpy as np

from ghist_secondary_contact.calibration import (
    apply_log_correction,
    competition_scores,
    fit_log_correction,
)


ESTIMATE_KEYS = (
    "population1_effective_size",
    "population2_effective_size",
    "generations_since_split",
    "generations_since_contact",
    "migration_rate",
)
TRUTH_COLUMNS = (1, 2, 3, 4, 5)
RIDGE_VALUES = (100.0, 30.0, 10.0, 3.0, 1.0, 0.3, 0.1)
SCALED_PARAMETER_KEYS = (
    "population1_size_ratio",
    "population2_size_ratio",
    "split_duration_2N",
    "contact_fraction",
    "migration_2N",
)
LOWER_BOUNDS = np.asarray((0.01, 0.01, 1e-5, 0.002, 1e-6))
UPPER_BOUNDS = np.asarray((100.0, 100.0, 10.0, 0.98, 100.0))


def score_summary(predictions: np.ndarray, truth: np.ndarray) -> dict:
    scores = competition_scores(predictions, truth)
    relative = predictions / truth - 1.0
    return {
        "replicates": len(scores),
        "mean_competition_rrmse": float(scores.mean()),
        "median_competition_rrmse": float(np.median(scores)),
        "competition_rrmse_interval_95": np.quantile(scores, (0.025, 0.975)).tolist(),
        "median_absolute_relative_error": dict(
            zip(ESTIMATE_KEYS, np.median(np.abs(relative), axis=0).tolist())
        ),
        "mean_relative_error": dict(zip(ESTIMATE_KEYS, relative.mean(axis=0).tolist())),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--design", type=Path, required=True)
    parser.add_argument("--inputs", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    design = json.loads(args.design.read_text())
    paths = sorted(Path(path) for path in glob.glob(args.inputs))
    records = {}
    for path in paths:
        record = json.loads(path.read_text())
        records[int(record["replicate"])] = record
    if len(records) != len(paths):
        raise ValueError("duplicate task indices in calibration results")
    expected = {int(row["task_index"]) for row in design["rows"]}
    if set(records) != expected:
        raise ValueError(
            f"result task indices differ from design: missing={sorted(expected - set(records))}, "
            f"extra={sorted(set(records) - expected)}"
        )
    boundary_tasks = []
    for index, record in records.items():
        scaled = np.asarray(
            [record["scaled_parameters"][key] for key in SCALED_PARAMETER_KEYS]
        )
        near_boundary = (scaled <= LOWER_BOUNDS * 1.001) | (
            scaled >= UPPER_BOUNDS / 1.001
        )
        if np.any(near_boundary):
            boundary_tasks.append(
                {
                    "task_index": index,
                    "parameters": [
                        key
                        for key, flagged in zip(SCALED_PARAMETER_KEYS, near_boundary)
                        if flagged
                    ],
                }
            )
    if boundary_tasks:
        raise ValueError(f"optimizer-boundary results require retry: {boundary_tasks}")

    rows = sorted(design["rows"], key=lambda row: int(row["task_index"]))
    estimates = np.asarray(
        [[records[index][key] for key in ESTIMATE_KEYS] for index in range(len(rows))],
        dtype=np.float64,
    )
    truth = np.asarray(
        [[row["truth"][column] for column in TRUTH_COLUMNS] for row in rows],
        dtype=np.float64,
    )
    partitions = np.asarray([row["partition"] for row in rows])
    masks = {name: partitions == name for name in ("calibration", "selection", "evaluation")}
    if any(not np.any(mask) for mask in masks.values()):
        raise ValueError("all three point-grouped partitions must be present")

    methods = ("identity", "median_ratio", *(f"ridge_{value:g}" for value in RIDGE_VALUES))
    selection_results = {}
    for method in methods:
        correction = fit_log_correction(
            method, estimates[masks["calibration"]], truth[masks["calibration"]]
        )
        predictions = apply_log_correction(correction, estimates[masks["selection"]])
        selection_results[method] = score_summary(predictions, truth[masks["selection"]])
    selected = min(
        methods,
        key=lambda method: (
            selection_results[method]["median_competition_rrmse"],
            selection_results[method]["mean_competition_rrmse"],
        ),
    )

    development_mask = masks["calibration"] | masks["selection"]
    final_correction = fit_log_correction(
        selected, estimates[development_mask], truth[development_mask]
    )
    evaluation_predictions = apply_log_correction(
        final_correction, estimates[masks["evaluation"]]
    )
    raw_evaluation = score_summary(estimates[masks["evaluation"]], truth[masks["evaluation"]])
    corrected_evaluation = score_summary(evaluation_predictions, truth[masks["evaluation"]])
    evaluation_improvement = (
        raw_evaluation["median_competition_rrmse"]
        - corrected_evaluation["median_competition_rrmse"]
    )
    raw_parameter_errors = raw_evaluation["median_absolute_relative_error"]
    corrected_parameter_errors = corrected_evaluation["median_absolute_relative_error"]
    parameter_error_ratios = {
        key: corrected_parameter_errors[key] / max(raw_parameter_errors[key], 1e-12)
        for key in ESTIMATE_KEYS
    }
    improved_parameters = sum(value < 1.0 for value in parameter_error_ratios.values())
    correction_accepted = (
        evaluation_improvement > 0
        and corrected_evaluation["mean_competition_rrmse"]
        < raw_evaluation["mean_competition_rrmse"]
        and improved_parameters >= 4
        and max(parameter_error_ratios.values()) <= 1.10
    )
    result = {
        "design": str(args.design),
        "inputs": args.inputs,
        "point_grouped_partitions": design["partition_point_counts"],
        "selection_rule": "lowest selection median competition RRMSE, then lowest mean",
        "selection_results": selection_results,
        "selected_method": selected,
        "untouched_evaluation": {
            "raw": raw_evaluation,
            "corrected": corrected_evaluation,
            "median_rrmse_improvement": evaluation_improvement,
            "parameter_error_ratios": parameter_error_ratios,
            "improved_parameters": improved_parameters,
            "acceptance_rule": "median and mean RRMSE improve, at least four parameter medians improve, and none worsens by more than 10%",
            "correction_accepted": correction_accepted,
        },
        "accepted_prediction_method": selected if correction_accepted else "identity",
        "final_correction": None
        if final_correction is None
        else {
            "center": final_correction.center.tolist(),
            "scale": final_correction.scale.tolist(),
            "coefficients": final_correction.coefficients.tolist(),
        },
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
