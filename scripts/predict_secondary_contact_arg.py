from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path

import numpy as np
import torch

from ghist_secondary_contact.features import inverse_transform_targets
from ghist_secondary_contact.model import SecondaryContactNet


OUTPUT_NAMES = (
    "mainland_population_size",
    "island_population_size",
    "generations_since_split",
    "generations_since_migration",
    "migration_rate",
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--checkpoints", required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--tiles-per-block", type=int, default=20)
    parser.add_argument("--batch-size", type=int, default=128)
    parser.add_argument(
        "--device",
        default="auto",
        choices=("auto", "cpu", "cuda"),
        help="Inference device; auto uses CUDA when it is available.",
    )
    args = parser.parse_args()
    device = torch.device(
        "cuda" if args.device == "auto" and torch.cuda.is_available() else args.device
    )
    if args.device == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA was requested but is not available")
    paths = [Path(path) for path in sorted(glob.glob(args.checkpoints))]
    if not paths:
        raise FileNotFoundError(f"no checkpoints match {args.checkpoints!r}")
    data = np.load(args.features)
    global_summary = np.concatenate(
        (
            data["singer_tile_mean"].reshape(len(data["singer_tile_mean"]), -1),
            data["singer_tile_posterior_sd"].reshape(
                len(data["singer_tile_posterior_sd"]), -1
            ),
            data["polegon_tile_mean"].reshape(len(data["polegon_tile_mean"]), -1),
        ),
        axis=1,
    )
    tiles = len(global_summary)
    if tiles % args.tiles_per_block:
        raise ValueError("tile count must be divisible by tiles per block")
    model_records = []
    model_block_predictions = []
    for path in paths:
        checkpoint = torch.load(path, map_location="cpu", weights_only=False)
        model = SecondaryContactNet(**checkpoint["model_configuration"])
        model.load_state_dict(checkpoint["model"])
        model.to(device).eval()
        normalization = checkpoint["normalization"]
        window_mean = np.asarray(normalization["window_mean"])
        window_scale = np.asarray(normalization["window_scale"])
        global_mean = np.asarray(normalization["global_mean"])
        global_scale = np.asarray(normalization["global_scale"])
        target_mean = np.asarray(normalization["target_mean"])
        target_scale = np.asarray(normalization["target_scale"])
        windows = (
            np.zeros((tiles, 1, checkpoint["model_configuration"]["window_features"]))
            - window_mean
        ) / window_scale
        spectra = np.zeros((tiles, 1, 87, 49), dtype=np.float32)
        globals_normalized = (global_summary - global_mean) / global_scale
        transformed = []
        log_variances = []
        with torch.no_grad():
            for start in range(0, tiles, args.batch_size):
                stop = min(start + args.batch_size, tiles)
                mean, log_variance = model(
                    torch.as_tensor(
                        windows[start:stop], dtype=torch.float32, device=device
                    ),
                    torch.as_tensor(
                        spectra[start:stop], dtype=torch.float32, device=device
                    ),
                    torch.as_tensor(
                        globals_normalized[start:stop],
                        dtype=torch.float32,
                        device=device,
                    ),
                )
                transformed.append(
                    mean.detach().cpu().numpy() * target_scale + target_mean
                )
                log_variances.append(log_variance.detach().cpu().numpy())
        transformed = np.concatenate(transformed)
        log_variances = np.concatenate(log_variances)
        block_prediction = transformed.reshape(
            -1, args.tiles_per_block, transformed.shape[1]
        ).mean(axis=1)
        model_center = np.median(block_prediction, axis=0)
        model_block_predictions.append(block_prediction)
        model_records.append(
            {
                "checkpoint": str(path),
                "best_epoch": checkpoint["epoch"],
                "selection_rrmse": checkpoint["selection_rrmse"],
                "prediction": dict(
                    zip(OUTPUT_NAMES, inverse_transform_targets(model_center).tolist())
                ),
                "median_tile_log_variance": np.median(log_variances, axis=0).tolist(),
            }
        )
    block_predictions = np.asarray(model_block_predictions)
    ensemble_blocks = block_predictions.mean(axis=0)
    center = np.median(ensemble_blocks, axis=0)
    physical_blocks = inverse_transform_targets(ensemble_blocks)
    result = {
        "policy": "local_prediction_only",
        "device": str(device),
        "tiles": tiles,
        "tiles_per_block": args.tiles_per_block,
        "blocks": len(ensemble_blocks),
        "models": model_records,
        "transformed_coordinate_ensemble": dict(
            zip(OUTPUT_NAMES, inverse_transform_targets(center).tolist())
        ),
        "block_percentile_05": dict(
            zip(OUTPUT_NAMES, np.percentile(physical_blocks, 5, axis=0).tolist())
        ),
        "block_percentile_95": dict(
            zip(OUTPUT_NAMES, np.percentile(physical_blocks, 95, axis=0).tolist())
        ),
        "model_spread_transformed": np.std(
            np.median(block_predictions, axis=1), axis=0
        ).tolist(),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
