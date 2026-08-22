from __future__ import annotations

import argparse
import json
from pathlib import Path

import joblib
import numpy as np
import torch

from ghist_bottleneck.infer import _robust_block_average
from ghist_bottleneck.model import BottleneckNet
from ghist_bottleneck.train import reference_features


def stratified_cell_split(
    cell_ids: np.ndarray,
    scenarios: np.ndarray,
    *,
    seed: int = 4421,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    calibration, selection, evaluation = [], [], []
    for scenario in np.unique(scenarios):
        group = np.asarray(cell_ids[scenarios == scenario]).copy()
        rng.shuffle(group)
        first = max(1, round(0.4 * len(group)))
        second = max(first + 1, round(0.7 * len(group)))
        second = min(second, len(group) - 1)
        calibration.extend(group[:first])
        selection.extend(group[first:second])
        evaluation.extend(group[second:])
    return tuple(np.asarray(values, dtype=np.int64) for values in (calibration, selection, evaluation))


def fit_affine(prediction: np.ndarray, target: np.ndarray, ridge: float = 1e-3) -> tuple[np.ndarray, np.ndarray]:
    design = np.column_stack((prediction, np.ones(len(prediction))))
    penalty = np.diag([ridge] * prediction.shape[1] + [0.0])
    coefficients = np.linalg.solve(design.T @ design + penalty, design.T @ target)
    return coefficients[:-1], coefficients[-1]


def apply_affine(prediction: np.ndarray, matrix: np.ndarray, intercept: np.ndarray) -> np.ndarray:
    return prediction @ matrix + intercept


def _network_predictions(
    windows: np.ndarray,
    spectra: np.ndarray,
    checkpoint_dir: Path,
    *,
    batch_size: int = 512,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    normalization = np.load(checkpoint_dir / "normalization.npz")
    normalized = ((windows - normalization["mean"]) / normalization["std"]).astype(np.float32)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    member_means, member_variances = [], []
    for path in sorted(checkpoint_dir.glob("model_*.pt")):
        state = torch.load(path, map_location=device)
        model = BottleneckNet(**state["model_config"]).to(device)
        model.load_state_dict(state["state_dict"])
        model.eval()
        means, variances = [], []
        with torch.inference_mode():
            for start in range(0, len(normalized), batch_size):
                stop = min(start + batch_size, len(normalized))
                mean, log_variance = model(
                    torch.from_numpy(normalized[start:stop]).to(device),
                    torch.from_numpy(spectra[start:stop].astype(np.float32)).to(device),
                )
                means.append(mean.cpu().numpy())
                variances.append(np.exp(log_variance.cpu().numpy()))
        member_means.append(
            np.concatenate(means) * normalization["target_std"] + normalization["target_mean"]
        )
        member_variances.append(
            np.concatenate(variances) * np.square(normalization["target_std"])
        )
    if not member_means:
        raise FileNotFoundError(f"no model checkpoints below {checkpoint_dir}")
    reference = joblib.load(checkpoint_dir / "reference.joblib")
    reference_prediction = reference["model"].predict(reference_features(windows, spectra))
    return (
        np.asarray(member_means),
        np.asarray(member_variances),
        reference_prediction,
        np.asarray(reference["weights"]),
    )


def _aggregate_cells(
    member_means: np.ndarray,
    member_variances: np.ndarray,
    reference_prediction: np.ndarray,
    blend: np.ndarray,
    cell_id: np.ndarray,
    n_variants: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    predictions, deep_variances = [], []
    for cell in np.unique(cell_id):
        mask = cell_id == cell
        deep_mean, deep_variance, _ = _robust_block_average(
            member_means[:, mask], member_variances[:, mask], n_variants[mask]
        )
        information = np.sqrt(n_variants[mask] / np.median(n_variants[mask]))
        information = np.clip(information, 0.25, 4.0)
        reference_mean = np.average(reference_prediction[mask], axis=0, weights=information)
        predictions.append(blend[0] * deep_mean + blend[1] * reference_mean)
        deep_variances.append(deep_variance)
    return np.asarray(predictions), np.asarray(deep_variances)


def _metrics(prediction: np.ndarray, target: np.ndarray, scale: np.ndarray) -> dict:
    residual = prediction - target
    return {
        "examples": len(target),
        "rmse_latent": np.sqrt(np.mean(np.square(residual), axis=0)).tolist(),
        "normalized_rmse": (np.sqrt(np.mean(np.square(residual), axis=0)) / scale).tolist(),
        "bias_latent": residual.mean(axis=0).tolist(),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", type=Path, required=True)
    parser.add_argument("--checkpoint-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    data = np.load(args.dataset)
    cell_ids = np.unique(data["cell_id"])
    first = np.asarray([np.flatnonzero(data["cell_id"] == cell)[0] for cell in cell_ids])
    target = data["targets"][first]
    scenarios = data["nuisance"][first, 0].astype(int)
    calibration_cells, selection_cells, evaluation_cells = stratified_cell_split(cell_ids, scenarios)
    index = {int(cell): position for position, cell in enumerate(cell_ids)}
    split_indices = {
        name: np.asarray([index[int(cell)] for cell in values], dtype=np.int64)
        for name, values in (
            ("calibration", calibration_cells),
            ("selection", selection_cells),
            ("evaluation", evaluation_cells),
        )
    }
    scale = target[split_indices["calibration"]].std(axis=0) + 1e-6
    branches = {}
    for name, spectra in (
        ("summary", np.array(data["sfs"], copy=True)),
        ("genealogy_proxy", np.array(data["sfs"], copy=True)),
    ):
        if name == "summary":
            spectra[:, -33:] = 0
        member_means, member_variances, reference_prediction, blend = _network_predictions(
            data["windows"], spectra, args.checkpoint_dir
        )
        raw, deep_variance = _aggregate_cells(
            member_means,
            member_variances,
            reference_prediction,
            blend,
            data["cell_id"],
            data["n_variants"],
        )
        calibration_index = split_indices["calibration"]
        matrix, intercept = fit_affine(raw[calibration_index], target[calibration_index])
        calibrated = apply_affine(raw, matrix, intercept)
        branches[name] = {
            "matrix": matrix.tolist(),
            "intercept": intercept.tolist(),
            "raw": raw,
            "calibrated": calibrated,
            "deep_variance": deep_variance,
            "selection": _metrics(
                calibrated[split_indices["selection"]], target[split_indices["selection"]], scale
            ),
            "evaluation": _metrics(
                calibrated[split_indices["evaluation"]], target[split_indices["evaluation"]], scale
            ),
        }
    # The simulated genealogy channel is a useful upper-bound diagnostic, but it
    # is not exchangeable with FastCXT's inferred hazards. Only the summary branch
    # is eligible until a matched translation-validation set exists.
    chosen = "summary"
    serializable = {
        "dataset": str(args.dataset),
        "checkpoint_dir": str(args.checkpoint_dir),
        "cells": len(cell_ids),
        "split_seed": 4421,
        "split_cells": {name: values.tolist() for name, values in zip(
            ("calibration", "selection", "evaluation"),
            (calibration_cells, selection_cells, evaluation_cells),
        )},
        "target_scale": scale.tolist(),
        "chosen_branch": chosen,
        "genealogy_proxy_eligible_for_observed_selection": False,
        "branches": {
            name: {key: value for key, value in branch.items() if key not in {"raw", "calibrated", "deep_variance"}}
            for name, branch in branches.items()
        },
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(serializable, indent=2) + "\n")
    np.savez_compressed(
        args.out.with_suffix(".npz"),
        cell_id=cell_ids,
        target_latent=target,
        scenario=scenarios,
        **{f"{name}_raw_latent": branch["raw"] for name, branch in branches.items()},
        **{f"{name}_calibrated_latent": branch["calibrated"] for name, branch in branches.items()},
    )
    print(json.dumps(serializable, indent=2))


if __name__ == "__main__":
    main()
