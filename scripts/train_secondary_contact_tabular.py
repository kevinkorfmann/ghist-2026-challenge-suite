from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path

import joblib
import numpy as np
from sklearn.decomposition import PCA
from sklearn.ensemble import ExtraTreesRegressor, HistGradientBoostingRegressor
from sklearn.linear_model import Ridge
from sklearn.multioutput import MultiOutputRegressor

from ghist_secondary_contact.features import (
    inverse_transform_targets,
    transform_targets,
)
from ghist_secondary_contact.contract import (
    SPECTRUM_CONTRACT,
    require_spectrum_contract,
)


def _load(pattern: str) -> dict[str, np.ndarray]:
    paths = [Path(path) for path in sorted(glob.glob(pattern))]
    if not paths:
        raise FileNotFoundError(f"no chunks match {pattern!r}")
    chunks = [np.load(path) for path in paths]
    for chunk in chunks:
        require_spectrum_contract(chunk)
    keys = ("partition", "design_domain", "truth", "point_id", "joint_sfs")
    result = {key: np.concatenate([chunk[key] for chunk in chunks]) for key in keys}
    result["spectrum_transform"] = np.concatenate(
        [
            np.repeat(str(chunk["spectrum_transform"]), len(chunk["truth"]))
            for chunk in chunks
        ]
    )
    result["n_variants"] = np.concatenate([chunk["n_variants"] for chunk in chunks])
    return result


def _point_rrmse(
    prediction: np.ndarray, truth: np.ndarray, point_id: np.ndarray
) -> float:
    estimates, targets = [], []
    for identifier in np.unique(point_id):
        selected = point_id == identifier
        estimates.append(prediction[selected].mean(axis=0))
        targets.append(truth[selected][0])
    relative = (np.asarray(estimates) - np.asarray(targets)) / np.asarray(targets)
    return float(np.sqrt(np.mean(np.sum(relative**2, axis=1))))


def _features(data: dict[str, np.ndarray]) -> np.ndarray:
    spectrum = data["joint_sfs"].reshape(len(data["truth"]), -1).astype(np.float32)
    count = np.log1p(data["n_variants"].astype(np.float64))[:, None].astype(np.float32)
    return np.concatenate((spectrum, count), axis=1)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--chunks", required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--domain-filter", required=True)
    parser.add_argument("--pca-components", type=int, default=96)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--seed", type=int, default=20260923)
    args = parser.parse_args()
    data = _load(args.chunks)
    domains = np.asarray(
        [int(value) for value in args.domain_filter.replace(":", ",").split(",")]
    )
    keep = np.isin(data["design_domain"], domains)
    data = {key: value[keep] for key, value in data.items()}
    transforms = np.unique(data["spectrum_transform"])
    if len(transforms) != 1:
        raise ValueError("Training chunks use different spectrum transforms")
    calibration = data["partition"] == 0
    selection = data["partition"] == 1
    evaluation = data["partition"] == 2
    if not all(np.any(values) for values in (calibration, selection, evaluation)):
        raise ValueError("all three grouped partitions must be present")
    raw_features = _features(data)
    components = min(args.pca_components, calibration.sum() - 1, raw_features.shape[1])
    pca = PCA(n_components=components, svd_solver="randomized", random_state=args.seed)
    pca.fit(raw_features[calibration])
    features = pca.transform(raw_features)
    transformed = transform_targets(data["truth"])
    target_mean = transformed[calibration].mean(axis=0)
    target_scale = transformed[calibration].std(axis=0).clip(1e-8)
    targets = (transformed - target_mean) / target_scale
    candidates = {
        "ridge_0.1": Ridge(alpha=0.1),
        "ridge_1": Ridge(alpha=1.0),
        "ridge_10": Ridge(alpha=10.0),
        "extra_trees_leaf_1": ExtraTreesRegressor(
            n_estimators=512,
            min_samples_leaf=1,
            max_features=0.8,
            n_jobs=args.workers,
            random_state=args.seed,
        ),
        "extra_trees_leaf_2": ExtraTreesRegressor(
            n_estimators=512,
            min_samples_leaf=2,
            max_features=0.8,
            n_jobs=args.workers,
            random_state=args.seed + 1,
        ),
        "extra_trees_leaf_4": ExtraTreesRegressor(
            n_estimators=512,
            min_samples_leaf=4,
            max_features=1.0,
            n_jobs=args.workers,
            random_state=args.seed + 2,
        ),
        "histogram_boosting": MultiOutputRegressor(
            HistGradientBoostingRegressor(
                learning_rate=0.05,
                max_iter=500,
                max_leaf_nodes=31,
                l2_regularization=0.1,
                early_stopping=True,
                random_state=args.seed,
            ),
            n_jobs=min(args.workers, 5),
        ),
    }
    selection_results = []
    best_name = None
    best_model = None
    best_score = np.inf
    for name, model in candidates.items():
        model.fit(features[calibration], targets[calibration])
        prediction = inverse_transform_targets(
            model.predict(features[selection]) * target_scale + target_mean
        )
        score = _point_rrmse(
            prediction,
            data["truth"][selection],
            data["point_id"][selection],
        )
        selection_results.append({"model": name, "selection_rrmse": score})
        print(json.dumps(selection_results[-1]), flush=True)
        if score < best_score:
            best_name, best_model, best_score = name, model, score
    assert best_name is not None and best_model is not None
    evaluation_prediction = inverse_transform_targets(
        best_model.predict(features[evaluation]) * target_scale + target_mean
    )
    evaluation_score = _point_rrmse(
        evaluation_prediction,
        data["truth"][evaluation],
        data["point_id"][evaluation],
    )
    report = {
        "selected_model": best_name,
        "selection_rrmse": best_score,
        "untouched_evaluation_rrmse": evaluation_score,
        "domain_filter": args.domain_filter,
        "pca_components": components,
        "pca_explained_variance": float(pca.explained_variance_ratio_.sum()),
        "candidate_selection_results": selection_results,
        "seed": args.seed,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(
        {
            "spectrum_contract": SPECTRUM_CONTRACT,
            "spectrum_transform": str(transforms[0]),
            "pca": pca,
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
