from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from sklearn.linear_model import LinearRegression, Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


def predict(model_name: str, train_x: np.ndarray, train_y: np.ndarray, test_x: np.ndarray) -> np.ndarray:
    if model_name == "identity":
        return np.exp(test_x)
    if model_name == "median_log_bias":
        correction = np.median(train_y - train_x, axis=0)
        return np.exp(test_x + correction)
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--recovery", type=Path, required=True)
    parser.add_argument("--observed", type=Path, action="append", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    recovery = json.loads(args.recovery.read_text())["results"]
    rows = list(recovery.values())
    raw = np.asarray(
        [
            [
                row["exponential"]["generations"],
                row["exponential"]["post_growth_scaling"],
            ]
            for row in rows
        ],
        dtype=np.float64,
    )
    truth = np.asarray(
        [
            [row["truth"]["generations"], row["truth"]["post_growth_scaling"]]
            for row in rows
        ],
        dtype=np.float64,
    )
    log_raw = np.log(raw)
    log_truth = np.log(truth)
    model_names = [
        "identity",
        "median_log_bias",
        "diagonal_affine",
        "joint_ridge_0.1",
        "joint_ridge_1.0",
        "joint_ridge_10.0",
    ]
    evaluations = {}
    for model_name in model_names:
        predictions = []
        for held_out in range(len(rows)):
            keep = np.arange(len(rows)) != held_out
            predictions.append(
                predict(model_name, log_raw[keep], log_truth[keep], log_raw[held_out])
            )
        predictions = np.asarray(predictions)
        relative_error = predictions / truth - 1
        per_truth_rrmse = np.sqrt(np.mean(relative_error**2, axis=1))
        evaluations[model_name] = {
            "mean_competition_rrmse": float(np.mean(per_truth_rrmse)),
            "median_competition_rrmse": float(np.median(per_truth_rrmse)),
            "maximum_competition_rrmse": float(np.max(per_truth_rrmse)),
            "generation_rmse": float(np.sqrt(np.mean(relative_error[:, 0] ** 2))),
            "scaling_rmse": float(np.sqrt(np.mean(relative_error[:, 1] ** 2))),
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

    observed_rows = []
    for path in args.observed:
        payload = json.loads(path.read_text())
        value = payload["results"]["noncoding_1kb_100kb"]
        estimate = np.asarray(
            [value["generations"], value["post_growth_scaling"]], dtype=np.float64
        )
        calibrated = predict(selected, log_raw, log_truth, np.log(estimate))
        observed_rows.append(
            {
                "projection": payload["projection"],
                "raw": {
                    "generations": float(estimate[0]),
                    "post_growth_scaling": float(estimate[1]),
                },
                "calibrated": {
                    "generations": float(calibrated[0]),
                    "post_growth_scaling": float(calibrated[1]),
                },
            }
        )
    calibrated_values = np.asarray(
        [
            [row["calibrated"]["generations"], row["calibrated"]["post_growth_scaling"]]
            for row in observed_rows
        ]
    )
    result = {
        "validation": "leave_one_truth_point_out",
        "truth_points": len(rows),
        "evaluations": evaluations,
        "selected_model": selected,
        "observed_by_projection": observed_rows,
        "calibrated_projection_median": {
            "generations": float(np.median(calibrated_values[:, 0])),
            "post_growth_scaling": float(np.median(calibrated_values[:, 1])),
        },
        "submission_count": 0,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
