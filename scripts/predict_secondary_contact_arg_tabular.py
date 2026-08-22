from __future__ import annotations

import argparse
import json
from pathlib import Path

import joblib
import numpy as np

from ghist_secondary_contact.features import inverse_transform_targets


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
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--tiles-per-block", type=int, default=20)
    args = parser.parse_args()

    data = np.load(args.features)
    features = np.concatenate(
        (
            data["singer_tile_mean"].reshape(len(data["singer_tile_mean"]), -1),
            data["singer_tile_posterior_sd"].reshape(
                len(data["singer_tile_posterior_sd"]), -1
            ),
            data["polegon_tile_mean"].reshape(len(data["polegon_tile_mean"]), -1),
        ),
        axis=1,
    ).astype(np.float32)
    features = np.nan_to_num(features)
    if len(features) % args.tiles_per_block:
        raise ValueError("tile count must be divisible by tiles per block")
    checkpoint = joblib.load(args.model)
    transformed = (
        checkpoint["model"].predict(features) * checkpoint["target_scale"]
        + checkpoint["target_mean"]
    )
    block_transformed = transformed.reshape(
        -1, args.tiles_per_block, transformed.shape[1]
    ).mean(axis=1)
    center = np.median(block_transformed, axis=0)
    block_prediction = inverse_transform_targets(block_transformed)
    result = {
        "policy": "local_prediction_only",
        "model": str(args.model),
        "selection_report": checkpoint["report"],
        "tiles": len(features),
        "tiles_per_block": args.tiles_per_block,
        "blocks": len(block_prediction),
        "prediction": dict(zip(OUTPUT_NAMES, inverse_transform_targets(center).tolist())),
        "block_percentile_05": dict(
            zip(OUTPUT_NAMES, np.percentile(block_prediction, 5, axis=0).tolist())
        ),
        "block_percentile_95": dict(
            zip(OUTPUT_NAMES, np.percentile(block_prediction, 95, axis=0).tolist())
        ),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
