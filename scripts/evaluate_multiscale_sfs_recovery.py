from __future__ import annotations

import argparse
import glob
import hashlib
import json
from pathlib import Path

import joblib
import numpy as np
from sklearn.ensemble import ExtraTreesRegressor, HistGradientBoostingRegressor
from sklearn.multioutput import MultiOutputRegressor
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import PolynomialFeatures, StandardScaler
from sklearn.linear_model import Ridge
from sklearn.neural_network import MLPRegressor

from run_multiscale_sfs_recovery import TOTAL_CELLS, TOTAL_POINTS


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def to_latent(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    fraction = np.clip(values[:, 1], 1e-6, 1 - 1e-6)
    return np.column_stack(
        (np.log(values[:, 0]), np.log(fraction / (1 - fraction)), np.log(values[:, 2]))
    )


def to_natural(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    return np.column_stack(
        (np.exp(values[:, 0]), 1 / (1 + np.exp(-values[:, 1])), np.exp(values[:, 2]))
    )


def metrics(
    prediction: np.ndarray,
    truth: np.ndarray,
    point_id: np.ndarray,
    scale: np.ndarray,
) -> dict:
    residual = prediction[:, :2] - truth[:, :2]
    prediction_natural = to_natural(prediction)
    truth_natural = to_natural(truth)
    relative = np.abs(prediction_natural[:, :2] / truth_natural[:, :2] - 1)
    point_bias = []
    for point in np.unique(point_id):
        mask = point_id == point
        point_bias.append(
            np.abs(prediction_natural[mask, :2].mean(axis=0) / truth_natural[mask, :2][0] - 1)
        )
    point_bias = np.asarray(point_bias)
    relative_rmse = np.sqrt(np.mean(np.square(relative), axis=0))
    return {
        "cells": len(prediction),
        "parameter_points": len(np.unique(point_id)),
        "normalized_rmse": (np.sqrt(np.mean(np.square(residual), axis=0)) / scale).tolist(),
        "median_absolute_relative_error": np.median(relative, axis=0).tolist(),
        "mean_absolute_relative_error": relative.mean(axis=0).tolist(),
        "relative_rmse": relative_rmse.tolist(),
        "rrmse": float(np.sqrt(np.mean(np.square(relative)))),
        "median_point_mean_absolute_relative_bias": np.median(point_bias, axis=0).tolist(),
        "maximum_point_mean_absolute_relative_bias": np.max(point_bias, axis=0).tolist(),
    }


def candidates(seed: int) -> dict:
    result = {"raw": None}
    for degree, alphas in ((1, (1e-4, 1e-2, 1.0)), (2, (1e-2, 0.1, 1.0, 10.0))):
        for alpha in alphas:
            result[f"ridge_d{degree}_a{alpha:g}"] = make_pipeline(
                PolynomialFeatures(degree=degree, include_bias=False),
                StandardScaler(),
                Ridge(alpha=alpha),
            )
    for leaf in (3, 6, 12):
        result[f"extra_trees_leaf{leaf}"] = ExtraTreesRegressor(
            n_estimators=500,
            min_samples_leaf=leaf,
            max_features=1.0,
            random_state=seed,
            n_jobs=1,
        )
    for regularization in (0.0, 0.1, 1.0):
        result[f"hist_gradient_l2_{regularization:g}"] = MultiOutputRegressor(
            HistGradientBoostingRegressor(
                max_iter=300,
                learning_rate=0.04,
                max_leaf_nodes=15,
                l2_regularization=regularization,
                random_state=seed,
            )
        )
    for alpha in (1e-4, 1e-2):
        result[f"mlp_alpha{alpha:g}"] = make_pipeline(
            StandardScaler(),
            MLPRegressor(
                hidden_layer_sizes=(64, 64),
                activation="tanh",
                alpha=alpha,
                learning_rate_init=2e-3,
                max_iter=1500,
                early_stopping=True,
                random_state=seed,
            ),
        )
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-glob", required=True)
    parser.add_argument("--observed", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--model-out", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20261201)
    args = parser.parse_args()
    paths = [Path(path) for path in sorted(glob.glob(args.input_glob))]
    if not paths:
        raise FileNotFoundError(args.input_glob)
    arrays = [np.load(path) for path in paths]
    fields = ("cell_id", "point_id", "replicate", "domain", "split", "truth", "estimate")
    merged = {field: np.concatenate([array[field] for array in arrays]) for field in fields}
    order = np.argsort(merged["cell_id"])
    merged = {field: values[order] for field, values in merged.items()}
    if len(merged["cell_id"]) != TOTAL_CELLS:
        raise ValueError(f"expected {TOTAL_CELLS} cells, found {len(merged['cell_id'])}")
    np.testing.assert_array_equal(merged["cell_id"], np.arange(TOTAL_CELLS))
    if len(np.unique(merged["point_id"])) != TOTAL_POINTS:
        raise ValueError("parameter point coverage is incomplete")
    for point in np.unique(merged["point_id"]):
        mask = merged["point_id"] == point
        if len(np.unique(merged["split"][mask])) != 1:
            raise ValueError("replicates from one parameter point crossed splits")

    truth = to_latent(merged["truth"])
    estimate = to_latent(merged["estimate"])
    calibration = merged["split"] == 0
    selection = merged["split"] == 1
    evaluation = merged["split"] == 2
    scale = truth[calibration, :2].std(axis=0) + 1e-6
    fitted = {}
    predictions = {}
    selection_metrics = {}
    raw_selection = {}
    for domain, name in ((0, "broad"), (1, "focused")):
        mask = selection & (merged["domain"] == domain)
        raw_selection[name] = metrics(
            estimate[mask], truth[mask], merged["point_id"][mask], scale
        )

    for name, model in candidates(args.seed).items():
        if model is None:
            prediction = np.array(estimate, copy=True)
        else:
            model.fit(estimate[calibration], truth[calibration])
            prediction = model.predict(estimate)
        fitted[name] = model
        predictions[name] = prediction
        per_domain = {}
        for domain, domain_name in ((0, "broad"), (1, "focused")):
            mask = selection & (merged["domain"] == domain)
            per_domain[domain_name] = metrics(
                prediction[mask], truth[mask], merged["point_id"][mask], scale
            )
        broad = per_domain["broad"]
        broad_raw = raw_selection["broad"]
        focused = per_domain["focused"]
        focused_raw = raw_selection["focused"]
        broad_improvement = broad["rrmse"] <= 0.98 * broad_raw["rrmse"]
        focused_guard = focused["rrmse"] <= 1.05 * focused_raw["rrmse"]
        per_domain["passes_selection_gate"] = bool(broad_improvement and focused_guard)
        per_domain["selection_score"] = float(
            (2 * broad["rrmse"] + focused["rrmse"]) / 3
        )
        selection_metrics[name] = per_domain

    eligible_selection = [
        name for name, value in selection_metrics.items() if value["passes_selection_gate"]
    ]
    selected = (
        min(eligible_selection, key=lambda name: selection_metrics[name]["selection_score"])
        if eligible_selection
        else "raw"
    )
    evaluation_metrics = {}
    for name in ("raw", selected):
        per_domain = {}
        for domain, domain_name in ((0, "broad"), (1, "focused")):
            mask = evaluation & (merged["domain"] == domain)
            per_domain[domain_name] = metrics(
                predictions[name][mask], truth[mask], merged["point_id"][mask], scale
            )
        evaluation_metrics[name] = per_domain
    raw_evaluation = evaluation_metrics["raw"]
    selected_evaluation = evaluation_metrics[selected]
    broad_improvement = (
        selected_evaluation["broad"]["rrmse"] <= 0.98 * raw_evaluation["broad"]["rrmse"]
    )
    focused_guard = (
        selected_evaluation["focused"]["rrmse"] <= 1.05 * raw_evaluation["focused"]["rrmse"]
    )
    eligible = bool(selected != "raw" and broad_improvement and focused_guard)

    observed = json.loads(args.observed.read_text())
    observed_natural = np.asarray(
        [[
            observed["generations"],
            observed["post_decline_fraction"],
            observed["ancestral_effective_size"],
        ]]
    )
    observed_latent = to_latent(observed_natural)
    observed_selected = (
        fitted[selected].predict(observed_latent) if selected != "raw" else observed_latent
    )
    observed_selected_natural = to_natural(observed_selected)[0]
    serializable = {
        "inputs": {path.name: sha256(path) for path in paths},
        "observed_sha256": sha256(args.observed),
        "evaluator_sha256": sha256(Path(__file__)),
        "cells": TOTAL_CELLS,
        "parameter_points": TOTAL_POINTS,
        "split_cells": {
            "calibration": int(calibration.sum()),
            "selection": int(selection.sum()),
            "evaluation": int(evaluation.sum()),
        },
        "selection_rule": "improve broad natural-scale RRMSE at least 2% while focused natural-scale RRMSE remains within 5% of raw; rank by two-thirds broad and one-third focused RRMSE",
        "final_rule": "repeat the broad RRMSE improvement and focused RRMSE guard on untouched parameter points",
        "selected_candidate": selected,
        "selection_metrics": selection_metrics,
        "evaluation_metrics": evaluation_metrics,
        "eligible": eligible,
        "chosen": selected if eligible else "raw",
        "observed": {
            "raw": observed_natural[0].tolist(),
            "selected_candidate": observed_selected_natural.tolist(),
            "chosen": (observed_selected_natural if eligible else observed_natural[0]).tolist(),
        },
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(serializable, indent=2) + "\n")
    args.model_out.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(
        {
            "selected_candidate": selected,
            "eligible": eligible,
            "model": fitted[selected],
            "input": "log generations, logit fraction, log ancestral size",
            "output": "same latent coordinates",
        },
        args.model_out,
    )
    print(json.dumps({key: value for key, value in serializable.items() if key != "selection_metrics"}, indent=2))


if __name__ == "__main__":
    main()
