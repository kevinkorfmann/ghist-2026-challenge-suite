from __future__ import annotations

from pathlib import Path
import json

import numpy as np
import torch
import yaml
import joblib

from .features import FeatureBundle, FeatureCollection
from .model import BottleneckNet
from .train import _device, reference_features


def apply_chromosome_calibration(
    latent: np.ndarray,
    calibration: str | Path | dict,
) -> tuple[np.ndarray, dict]:
    """Apply a calibration frozen by the chromosome-level validation gate."""
    record = (
        json.loads(Path(calibration).read_text())
        if not isinstance(calibration, dict)
        else calibration
    )
    branch_name = record["chosen_branch"]
    branch = record["branches"][branch_name]
    matrix = np.asarray(branch["matrix"], dtype=float)
    intercept = np.asarray(branch["intercept"], dtype=float)
    latent = np.asarray(latent, dtype=float)
    if latent.shape != (2,) or matrix.shape != (2, 2) or intercept.shape != (2,):
        raise ValueError("chromosome calibration must map one two-parameter latent vector")
    return latent @ matrix + intercept, {
        "source": str(calibration) if not isinstance(calibration, dict) else "in-memory",
        "branch": branch_name,
        "matrix": matrix.tolist(),
        "intercept": intercept.tolist(),
    }


def predict(observed: str | Path, checkpoint_dir: str | Path, out: str | Path) -> dict:
    b = FeatureBundle.load(observed)
    checkpoint_dir = Path(checkpoint_dir)
    norm = np.load(checkpoint_dir / "normalization.npz")
    xw = ((b.windows[None] - norm["mean"]) / norm["std"]).astype(np.float32)
    xs = b.sfs[None].astype(np.float32)
    device = _device()
    preds, aleatoric = [], []
    for path in sorted(checkpoint_dir.glob("model_*.pt")):
        state = torch.load(path, map_location=device)
        model_config = state.get(
            "model_config",
            {"window_features": b.windows.shape[1], "sfs_bins": b.sfs.shape[0]},
        )
        model = BottleneckNet(**model_config).to(device)
        model.load_state_dict(state["state_dict"])
        model.eval()
        with torch.no_grad():
            mean, log_var = model(torch.from_numpy(xw).to(device), torch.from_numpy(xs).to(device))
        preds.append(mean.cpu().numpy()[0])
        aleatoric.append(np.exp(log_var.cpu().numpy()[0]))
    if not preds:
        raise FileNotFoundError(f"No model_*.pt checkpoints in {checkpoint_dir}")
    z = np.asarray(preds)
    deep_mean = z.mean(0) * norm["target_std"][0] + norm["target_mean"][0]
    total_var = (np.asarray(aleatoric).mean(0) + z.var(0)) * np.square(norm["target_std"][0])
    mean = deep_mean
    component_predictions = {"deep_ensemble_latent": deep_mean.tolist()}
    reference_path = checkpoint_dir / "reference.joblib"
    if reference_path.exists():
        reference = joblib.load(reference_path)
        raw_w = b.windows[None].astype(np.float32)
        ref_mean = reference["model"].predict(reference_features(raw_w, xs))[0]
        weights = np.asarray(reference["weights"])
        mean = weights[0] * deep_mean + weights[1] * ref_mean
        component_predictions["reference_latent"] = ref_mean.tolist()
        component_predictions["blend_weights_deep_reference"] = weights.T.tolist()
    generations = float(np.exp(mean[0]))
    fraction = float(1 / (1 + np.exp(-mean[1])))
    result = {
        "generations": int(round(generations)),
        "post_decline_fraction": round(fraction, 6),
        "diagnostics": {
            "ensemble_members": len(preds),
            "latent_sd": [float(x) for x in np.sqrt(total_var)],
            "n_variants": b.n_variants,
            "n_haplotypes": b.n_haplotypes,
            "components": component_predictions,
        },
    }
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(yaml.safe_dump({k: result[k] for k in ("generations", "post_decline_fraction")}, sort_keys=False))
    out.with_suffix(out.suffix + ".json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def _robust_block_average(
    values: np.ndarray,
    variances: np.ndarray,
    variant_counts: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Aggregate matched genomic blocks while limiting local outlier leverage."""
    model_means, model_variances, block_scales = [], [], []
    valid = variant_counts > 0
    if not np.any(valid):
        raise ValueError("all observed blocks are empty")
    information = np.sqrt(variant_counts[valid] / np.median(variant_counts[valid]))
    information = np.clip(information, 0.25, 4.0)
    for member_values, member_variances in zip(values[:, valid], variances[:, valid]):
        center = np.median(member_values, axis=0)
        absolute = np.abs(member_values - center)
        scale = 1.4826 * np.median(absolute, axis=0) + 1e-6
        huber = np.minimum(1.0, 2.5 * scale / np.maximum(absolute, 1e-6))
        precision = 1 / np.maximum(member_variances, 1e-6)
        lower, upper = np.quantile(precision, [0.1, 0.9], axis=0)
        precision = np.clip(precision, lower, upper)
        weights = huber * precision * information[:, None]
        mean = np.sum(weights * member_values, axis=0) / np.sum(weights, axis=0)
        effective = np.square(np.sum(weights, axis=0)) / np.sum(np.square(weights), axis=0)
        between = np.sum(weights * np.square(member_values - mean), axis=0) / np.sum(
            weights, axis=0
        )
        variance = between / np.maximum(effective, 1.0) + 1 / np.sum(precision, axis=0)
        model_means.append(mean)
        model_variances.append(variance)
        block_scales.append(scale)
    model_means = np.asarray(model_means)
    model_variances = np.asarray(model_variances)
    ensemble_mean = model_means.mean(axis=0)
    total_variance = model_variances.mean(axis=0) + model_means.var(axis=0)
    return ensemble_mean, total_variance, np.asarray(block_scales)


def predict_blocks(
    observed: str | Path,
    checkpoint_dir: str | Path,
    out: str | Path,
    calibration: str | Path | dict | None = None,
) -> dict:
    collection = FeatureCollection.load(observed)
    checkpoint_dir = Path(checkpoint_dir)
    norm = np.load(checkpoint_dir / "normalization.npz")
    xw = ((collection.windows - norm["mean"]) / norm["std"]).astype(np.float32)
    xs = collection.sfs.astype(np.float32)
    device = _device()
    member_predictions, member_variances = [], []
    for path in sorted(checkpoint_dir.glob("model_*.pt")):
        state = torch.load(path, map_location=device)
        model_config = state.get(
            "model_config",
            {
                "window_features": collection.windows.shape[2],
                "sfs_bins": collection.sfs.shape[1],
            },
        )
        model = BottleneckNet(**model_config).to(device)
        model.load_state_dict(state["state_dict"])
        model.eval()
        with torch.no_grad():
            mean, log_var = model(
                torch.from_numpy(xw).to(device), torch.from_numpy(xs).to(device)
            )
        member_predictions.append(
            mean.cpu().numpy() * norm["target_std"] + norm["target_mean"]
        )
        member_variances.append(
            np.exp(log_var.cpu().numpy()) * np.square(norm["target_std"])
        )
    if not member_predictions:
        raise FileNotFoundError(f"No model_*.pt checkpoints in {checkpoint_dir}")
    deep_mean, total_variance, block_scales = _robust_block_average(
        np.asarray(member_predictions),
        np.asarray(member_variances),
        collection.n_variants,
    )
    mean = deep_mean
    components = {"deep_block_ensemble_latent": deep_mean.tolist()}
    reference_path = checkpoint_dir / "reference.joblib"
    reference_blocks = None
    if reference_path.exists():
        reference = joblib.load(reference_path)
        reference_blocks = reference["model"].predict(
            reference_features(collection.windows.astype(np.float32), xs)
        )
        valid = collection.n_variants > 0
        info = np.sqrt(collection.n_variants[valid] / np.median(collection.n_variants[valid]))
        info = np.clip(info, 0.25, 4.0)
        reference_mean = np.average(reference_blocks[valid], axis=0, weights=info)
        blend = np.asarray(reference["weights"])
        mean = blend[0] * deep_mean + blend[1] * reference_mean
        components["reference_block_latent"] = reference_mean.tolist()
        components["blend_weights_deep_reference"] = blend.T.tolist()
    calibration_record = None
    components["pre_calibration_latent"] = np.asarray(mean).tolist()
    if calibration is not None:
        mean, calibration_record = apply_chromosome_calibration(mean, calibration)
        components["calibrated_latent"] = np.asarray(mean).tolist()
    generations = float(np.exp(mean[0]))
    fraction = float(1 / (1 + np.exp(-mean[1])))
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(
        yaml.safe_dump(
            {
                "generations": int(round(generations)),
                "post_decline_fraction": round(fraction, 6),
            },
            sort_keys=False,
        )
    )
    diagnostics_path = out.with_suffix(out.suffix + ".blocks.npz")
    np.savez_compressed(
        diagnostics_path,
        deep_member_block_latent=np.asarray(member_predictions),
        deep_member_block_variance_latent=np.asarray(member_variances),
        reference_block_latent=reference_blocks
        if reference_blocks is not None
        else np.empty((0, 2)),
        block_starts=collection.block_starts,
        block_ends=collection.block_ends,
        n_variants=collection.n_variants,
    )
    result = {
        "generations": int(round(generations)),
        "post_decline_fraction": round(fraction, 6),
        "diagnostics": {
            "blocks_total": int(len(collection.n_variants)),
            "blocks_used": int(np.sum(collection.n_variants > 0)),
            "ensemble_members": len(member_predictions),
            "latent_sd": np.sqrt(total_variance).tolist(),
            "median_member_block_scale": np.median(block_scales, axis=0).tolist(),
            "components": components,
            "chromosome_calibration": calibration_record,
            "block_artifact": str(diagnostics_path),
        },
    }
    out.with_suffix(out.suffix + ".json").write_text(json.dumps(result, indent=2) + "\n")
    return result
