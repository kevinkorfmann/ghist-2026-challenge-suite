from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from sklearn.linear_model import LinearRegression, Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from ghist_growth_bgs.simulation import require_physical_scaling


PARAMETERS = (
    "AFR_population_size",
    "EUR_population_size",
    "generations_since_split",
    "migration_rate",
)


def predict(
    model_name: str,
    train_x: np.ndarray,
    train_y: np.ndarray,
    test_x: np.ndarray,
) -> np.ndarray:
    if model_name == "identity":
        return np.exp(test_x)
    if model_name == "median_log_bias":
        return np.exp(test_x + np.median(train_y - train_x, axis=0))
    if model_name == "diagonal_affine":
        values = []
        for column in range(train_x.shape[1]):
            model = LinearRegression().fit(train_x[:, [column]], train_y[:, column])
            values.append(model.predict(test_x[[column]].reshape(1, 1))[0])
        return np.exp(values)
    if model_name.startswith("joint_ridge_"):
        alpha = float(model_name.rsplit("_", 1)[1])
        model = make_pipeline(StandardScaler(), Ridge(alpha=alpha))
        model.fit(train_x, train_y)
        return np.exp(model.predict(test_x.reshape(1, -1))[0])
    raise ValueError(f"unknown calibration model: {model_name}")


def vector(row: dict) -> np.ndarray:
    return np.asarray([row[key] for key in PARAMETERS], dtype=np.float64)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--recovery", type=Path, action="append", required=True)
    parser.add_argument("--observed", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    recoveries = [json.loads(path.read_text()) for path in args.recovery]
    scaling_scheme = require_physical_scaling(recoveries)
    raw = np.asarray([vector(row["fit"]) for row in recoveries])
    truth = np.asarray([vector(row["truth"]) for row in recoveries])
    log_raw, log_truth = np.log(raw), np.log(truth)
    model_names = (
        "identity",
        "median_log_bias",
        "diagonal_affine",
        "joint_ridge_0.1",
        "joint_ridge_1.0",
        "joint_ridge_10.0",
    )
    evaluations = {}
    for model_name in model_names:
        predictions = []
        for held_out in range(len(recoveries)):
            keep = np.arange(len(recoveries)) != held_out
            predictions.append(
                predict(model_name, log_raw[keep], log_truth[keep], log_raw[held_out])
            )
        predictions = np.asarray(predictions)
        relative_error = predictions / truth - 1
        per_truth = np.sqrt(np.mean(relative_error**2, axis=1))
        evaluations[model_name] = {
            "mean_competition_rrmse": float(np.mean(per_truth)),
            "median_competition_rrmse": float(np.median(per_truth)),
            "maximum_competition_rrmse": float(np.max(per_truth)),
            "parameter_rmse": {
                key: float(np.sqrt(np.mean(relative_error[:, index] ** 2)))
                for index, key in enumerate(PARAMETERS)
            },
            "predictions": predictions.tolist(),
            "truth": truth.tolist(),
        }
    selected = min(
        model_names,
        key=lambda name: (
            evaluations[name]["mean_competition_rrmse"],
            evaluations[name]["maximum_competition_rrmse"],
        ),
    )
    observed_payload = json.loads(args.observed.read_text())
    observed = vector(observed_payload["projection_median"])
    calibrated = predict(selected, log_raw, log_truth, np.log(observed))
    payload = {
        "simulation_scaling_scheme": scaling_scheme,
        "validation": "leave_one_truth_point_out",
        "truth_points": len(recoveries),
        "parameters": list(PARAMETERS),
        "evaluations": evaluations,
        "selected_model": selected,
        "observed_raw": dict(zip(PARAMETERS, observed.tolist())),
        "observed_calibrated": dict(zip(PARAMETERS, calibrated.tolist())),
        "submission_count": 0,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
