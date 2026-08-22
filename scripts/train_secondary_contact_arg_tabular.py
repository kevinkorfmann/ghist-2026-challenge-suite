from __future__ import annotations

import argparse
import json
from pathlib import Path

import joblib
import numpy as np
from sklearn.ensemble import ExtraTreesRegressor, HistGradientBoostingRegressor
from sklearn.linear_model import Ridge
from sklearn.multioutput import MultiOutputRegressor
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from ghist_secondary_contact.features import inverse_transform_targets, transform_targets


def point_rrmse(
    prediction: np.ndarray, truth: np.ndarray, point_id: np.ndarray
) -> float:
    estimates, targets = [], []
    for identifier in np.unique(point_id):
        selected = point_id == identifier
        estimates.append(prediction[selected].mean(axis=0))
        targets.append(truth[selected][0])
    relative = (np.asarray(estimates) - np.asarray(targets)) / np.asarray(targets)
    return float(np.sqrt(np.mean(np.sum(relative**2, axis=1))))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--domain-filter", required=True)
    parser.add_argument(
        "--feature-key",
        choices=("global_summary", "truth_hazards"),
        default="global_summary",
    )
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--seed", type=int, default=20260929)
    args = parser.parse_args()

    source = np.load(args.dataset)
    domains = np.asarray(
        [int(value) for value in args.domain_filter.replace(":", ",").split(",")]
    )
    keep = np.isin(source["design_domain"], domains)
    data = {
        key: source[key][keep]
        for key in ("partition", "truth", "point_id", args.feature_key)
    }
    calibration = data["partition"] == 0
    selection = data["partition"] == 1
    evaluation = data["partition"] == 2
    if not all(np.any(values) for values in (calibration, selection, evaluation)):
        raise ValueError("all three grouped partitions must be present")

    features = np.nan_to_num(
        data[args.feature_key].reshape(len(data["truth"]), -1).astype(np.float32)
    )
    transformed = transform_targets(data["truth"])
    target_mean = transformed[calibration].mean(axis=0)
    target_scale = transformed[calibration].std(axis=0).clip(1e-8)
    targets = (transformed - target_mean) / target_scale
    candidates = {
        "ridge_0.1": make_pipeline(StandardScaler(), Ridge(alpha=0.1)),
        "ridge_1": make_pipeline(StandardScaler(), Ridge(alpha=1.0)),
        "ridge_10": make_pipeline(StandardScaler(), Ridge(alpha=10.0)),
        "extra_trees_leaf_2": ExtraTreesRegressor(
            n_estimators=1024,
            min_samples_leaf=2,
            max_features=0.8,
            n_jobs=args.workers,
            random_state=args.seed,
        ),
        "extra_trees_leaf_4": ExtraTreesRegressor(
            n_estimators=1024,
            min_samples_leaf=4,
            max_features=1.0,
            n_jobs=args.workers,
            random_state=args.seed + 1,
        ),
        "histogram_boosting": MultiOutputRegressor(
            HistGradientBoostingRegressor(
                learning_rate=0.04,
                max_iter=700,
                max_leaf_nodes=31,
                l2_regularization=0.2,
                early_stopping=True,
                random_state=args.seed,
            ),
            n_jobs=min(args.workers, 5),
        ),
    }
    selection_results = []
    best_name, best_model, best_score = None, None, np.inf
    for name, model in candidates.items():
        model.fit(features[calibration], targets[calibration])
        prediction = inverse_transform_targets(
            model.predict(features[selection]) * target_scale + target_mean
        )
        score = point_rrmse(
            prediction, data["truth"][selection], data["point_id"][selection]
        )
        record = {"model": name, "selection_rrmse": score}
        selection_results.append(record)
        print(json.dumps(record), flush=True)
        if score < best_score:
            best_name, best_model, best_score = name, model, score
    assert best_name is not None and best_model is not None
    evaluation_prediction = inverse_transform_targets(
        best_model.predict(features[evaluation]) * target_scale + target_mean
    )
    evaluation_score = point_rrmse(
        evaluation_prediction,
        data["truth"][evaluation],
        data["point_id"][evaluation],
    )
    report = {
        "selected_model": best_name,
        "selection_rrmse": best_score,
        "untouched_evaluation_rrmse": evaluation_score,
        "domain_filter": args.domain_filter,
        "feature_key": args.feature_key,
        "features": int(features.shape[1]),
        "calibration_rows": int(calibration.sum()),
        "selection_rows": int(selection.sum()),
        "evaluation_rows": int(evaluation.sum()),
        "candidate_selection_results": selection_results,
        "seed": args.seed,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(
        {
            "model": best_model,
            "target_mean": target_mean,
            "target_scale": target_scale,
            "report": report,
        },
        args.out,
    )
    args.out.with_suffix(".json").write_text(json.dumps(report, indent=2) + "\n")
    np.savez_compressed(
        args.out.with_suffix(".evaluation.npz"),
        indices=np.flatnonzero(evaluation),
        prediction=evaluation_prediction,
        truth=data["truth"][evaluation],
        point_id=data["point_id"][evaluation],
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
