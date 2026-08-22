from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--tiles", type=int, default=1000)
    parser.add_argument("--tiles-per-block", type=int, default=20)
    args = parser.parse_args()
    if args.tiles % args.tiles_per_block:
        raise ValueError("tiles must be divisible by tiles per block")
    singer_mean, singer_sd, polegon = [], [], []
    polegon_interval_counts = []
    missing = []
    for tile in range(args.tiles):
        root = args.root / f"tile_{tile:04d}"
        singer_path = root / "coalescence_hazards.npz"
        polegon_path = root / "polegon_coalescence_hazards.npz"
        if not singer_path.exists() or not polegon_path.exists():
            missing.append(tile)
            continue
        singer = np.load(singer_path)
        singer_mean.append(singer["mean"])
        singer_sd.append(singer["standard_deviation"])
        polegon_data = np.load(polegon_path)
        polegon.append(polegon_data["mean"])
        if "total_interval_counts" in polegon_data:
            polegon_interval_counts.append(polegon_data["total_interval_counts"])
    if missing:
        raise FileNotFoundError(f"missing tile outputs: {missing[:20]}")
    singer_mean = np.asarray(singer_mean)
    singer_sd = np.asarray(singer_sd)
    polegon = np.asarray(polegon)
    if polegon_interval_counts and len(polegon_interval_counts) != args.tiles:
        raise ValueError("interval counts are present for only a subset of tiles")
    blocks = args.tiles // args.tiles_per_block
    singer_blocks = singer_mean.reshape(blocks, args.tiles_per_block, 3, 32).mean(axis=1)
    polegon_blocks = polegon.reshape(blocks, args.tiles_per_block, 3, 32).mean(axis=1)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    arrays = dict(
        singer_tile_mean=singer_mean,
        singer_tile_posterior_sd=singer_sd,
        polegon_tile_mean=polegon,
        singer_block_mean=singer_blocks,
        polegon_block_mean=polegon_blocks,
        singer_global_mean=singer_mean.mean(axis=0),
        polegon_global_mean=polegon.mean(axis=0),
    )
    if polegon_interval_counts:
        count_array = np.asarray(polegon_interval_counts)
        arrays.update(
            polegon_tile_interval_counts=count_array,
            polegon_block_interval_counts=count_array.reshape(
                blocks, args.tiles_per_block, 3, 33
            ).sum(axis=1),
            polegon_global_interval_counts=count_array.sum(axis=0),
        )
    np.savez_compressed(args.out, **arrays)
    report = {
        "tiles": args.tiles,
        "tiles_per_block": args.tiles_per_block,
        "blocks": blocks,
        "singer_posterior_draws_per_tile": 4,
        "mean_absolute_redating_difference": float(np.mean(np.abs(singer_mean - polegon))),
        "pair_class_redating_difference": np.mean(
            np.abs(singer_mean - polegon), axis=(0, 2)
        ).tolist(),
        "median_singer_posterior_sd": np.median(singer_sd, axis=(0, 2)).tolist(),
        "polegon_interval_counts_available": bool(polegon_interval_counts),
    }
    args.out.with_suffix(args.out.suffix + ".json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
