from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

try:
    from scripts.evaluate_bottleneck_chromosomes import _metrics
except ModuleNotFoundError:
    from evaluate_bottleneck_chromosomes import _metrics


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def natural(latent: np.ndarray) -> np.ndarray:
    latent = np.asarray(latent, dtype=float)
    return np.column_stack((np.exp(latent[:, 0]), 1 / (1 + np.exp(-latent[:, 1]))))


def median_errors(prediction: np.ndarray, target: np.ndarray) -> np.ndarray:
    prediction_natural = natural(prediction)
    target_natural = natural(target)
    return np.array(
        [
            np.median(np.abs(prediction_natural[:, 0] / target_natural[:, 0] - 1)),
            np.median(np.abs(prediction_natural[:, 1] - target_natural[:, 1])),
        ]
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--moments", type=Path, required=True)
    parser.add_argument("--learned", type=Path, required=True)
    parser.add_argument("--split-record", type=Path, required=True)
    parser.add_argument("--observed-moments", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    moments = np.load(args.moments)
    learned = np.load(args.learned)
    split_record = json.loads(args.split_record.read_text())
    moments_index = {int(cell): index for index, cell in enumerate(moments["cell_id"])}
    learned_index = {int(cell): index for index, cell in enumerate(learned["cell_id"])}
    common = sorted(set(moments_index) & set(learned_index))
    target = np.asarray([moments["target_latent"][moments_index[cell]] for cell in common])
    mechanistic = np.asarray([moments["raw_latent"][moments_index[cell]] for cell in common])
    learned_prediction = np.asarray(
        [learned["summary_calibrated_latent"][learned_index[cell]] for cell in common]
    )
    index = {cell: position for position, cell in enumerate(common)}
    split_indices = {
        name: np.asarray([index[int(cell)] for cell in cells if int(cell) in index])
        for name, cells in split_record["split_cells"].items()
    }
    calibration = split_indices["calibration"]
    selection = split_indices["selection"]
    evaluation = split_indices["evaluation"]
    scale = target[calibration].std(axis=0) + 1e-6
    selection_raw_rmse = np.sqrt(
        np.mean(np.square(mechanistic[selection] - target[selection]), axis=0)
    ) / scale
    selection_raw_median = median_errors(mechanistic[selection], target[selection])
    weights = np.linspace(0, 1, 21)
    chosen_weights = []
    selection_candidates = []
    for coordinate in range(2):
        candidates = []
        for weight in weights:
            prediction = (
                weight * mechanistic[:, coordinate]
                + (1 - weight) * learned_prediction[:, coordinate]
            )
            rmse = float(
                np.sqrt(np.mean(np.square(prediction[selection] - target[selection, coordinate])))
                / scale[coordinate]
            )
            combined = np.array(mechanistic, copy=True)
            combined[:, coordinate] = prediction
            median = float(median_errors(combined[selection], target[selection])[coordinate])
            candidates.append(
                {
                    "mechanistic_weight": float(weight),
                    "normalized_rmse": rmse,
                    "natural_median_error": median,
                    "passes_natural_guard": bool(median <= 1.05 * selection_raw_median[coordinate]),
                }
            )
        eligible = [candidate for candidate in candidates if candidate["passes_natural_guard"]]
        winner = min(eligible, key=lambda candidate: candidate["normalized_rmse"])
        if winner["normalized_rmse"] > 0.98 * selection_raw_rmse[coordinate]:
            winner = next(candidate for candidate in candidates if candidate["mechanistic_weight"] == 1)
        chosen_weights.append(winner["mechanistic_weight"])
        selection_candidates.append(candidates)

    chosen_weights = np.asarray(chosen_weights)
    ensemble = chosen_weights * mechanistic + (1 - chosen_weights) * learned_prediction
    raw_metrics = {
        "selection": _metrics(mechanistic[selection], target[selection], scale, bootstrap_seed=8611),
        "evaluation": _metrics(mechanistic[evaluation], target[evaluation], scale, bootstrap_seed=8612),
    }
    ensemble_metrics = {
        "selection": _metrics(ensemble[selection], target[selection], scale, bootstrap_seed=8621),
        "evaluation": _metrics(ensemble[evaluation], target[evaluation], scale, bootstrap_seed=8622),
    }
    raw_evaluation = raw_metrics["evaluation"]
    ensemble_evaluation = ensemble_metrics["evaluation"]
    improves_rmse = np.all(
        np.asarray(ensemble_evaluation["normalized_rmse"])
        <= 0.98 * np.asarray(raw_evaluation["normalized_rmse"])
    )
    natural_guard = (
        ensemble_evaluation["generation_median_absolute_relative_error"]
        <= 1.05 * raw_evaluation["generation_median_absolute_relative_error"]
        and ensemble_evaluation["fraction_median_absolute_error"]
        <= 1.05 * raw_evaluation["fraction_median_absolute_error"]
    )
    eligible = bool(improves_rmse and natural_guard)

    observed_moments = json.loads(args.observed_moments.read_text())
    observed_mechanistic = np.array(
        [
            np.log(observed_moments["generations"]),
            np.log(
                observed_moments["post_decline_fraction"]
                / (1 - observed_moments["post_decline_fraction"])
            ),
        ]
    )
    observed_learned_record = split_record["observed_summary"]
    observed_learned = np.array(
        [
            np.log(observed_learned_record["generations"]),
            np.log(
                observed_learned_record["post_decline_fraction"]
                / (1 - observed_learned_record["post_decline_fraction"])
            ),
        ]
    )
    observed_ensemble = chosen_weights * observed_mechanistic + (1 - chosen_weights) * observed_learned
    observed_ensemble_natural = natural(observed_ensemble[None])[0]
    serializable = {
        "moments_sha256": sha256(args.moments),
        "learned_sha256": sha256(args.learned),
        "split_record_sha256": sha256(args.split_record),
        "observed_moments_sha256": sha256(args.observed_moments),
        "evaluator_sha256": sha256(Path(__file__)),
        "common_nominal_cells": len(common),
        "split_sizes": {name: len(values) for name, values in split_indices.items()},
        "selection_rule": "per-coordinate 0.05 weight grid; minimize selection normalized RMSE subject to natural median within 5% of raw; require at least 2% RMSE improvement",
        "final_gate": "both evaluation normalized RMSE values improve at least 2% and both natural medians remain within 5% of raw",
        "mechanistic_weights": chosen_weights.tolist(),
        "selection_candidates": selection_candidates,
        "raw_metrics": raw_metrics,
        "ensemble_metrics": ensemble_metrics,
        "eligible": eligible,
        "chosen": "ensemble" if eligible else "mechanistic",
        "observed": {
            "mechanistic_generations": float(np.exp(observed_mechanistic[0])),
            "mechanistic_post_decline_fraction": float(natural(observed_mechanistic[None])[0, 1]),
            "learned_generations": float(np.exp(observed_learned[0])),
            "learned_post_decline_fraction": float(natural(observed_learned[None])[0, 1]),
            "ensemble_generations": float(observed_ensemble_natural[0]),
            "ensemble_post_decline_fraction": float(observed_ensemble_natural[1]),
        },
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(serializable, indent=2) + "\n")
    print(json.dumps({key: value for key, value in serializable.items() if key != "selection_candidates"}, indent=2))


if __name__ == "__main__":
    main()
