from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from ghist_bottleneck.fastcxt_dataset import choose_pairs


def select_checkpoint(root: Path) -> tuple[Path, dict]:
    metric_paths = sorted(root.glob("**/metrics.csv"))
    if len(metric_paths) != 1:
        raise ValueError(f"expected one metrics.csv below {root}, found {len(metric_paths)}")
    rows = list(csv.DictReader(metric_paths[0].open()))
    validation = [row for row in rows if row.get("val_loss") not in (None, "")]
    if not validation:
        raise ValueError("training log has no validation losses")
    best = min(validation, key=lambda row: float(row["val_loss"]))
    epoch = int(float(best["epoch"]))
    checkpoints = sorted(root.glob(f"**/epoch={epoch}-step=*.ckpt"))
    if len(checkpoints) != 1:
        raise ValueError(f"expected one checkpoint for epoch {epoch}, found {len(checkpoints)}")
    return checkpoints[0], best


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix-dir", type=Path, required=True)
    parser.add_argument("--checkpoint-root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--pairs", type=int, default=128)
    parser.add_argument("--pair-seed", type=int, default=2_026_082_2)
    parser.add_argument("--mutation-rate", type=float, default=3.42e-8)
    parser.add_argument("--sequence-length", type=int, default=100_000_000)
    parser.add_argument("--block-length", type=int, default=1_000_000)
    args = parser.parse_args()

    from fastcxt.train import LitFastCxt
    from fastcxt.translate import translate_from_genotype_matrix

    checkpoint, validation = select_checkpoint(args.checkpoint_root)
    fitted = LitFastCxt.load_from_checkpoint(checkpoint, map_location="cpu")
    model = fitted.model
    genotypes = np.load(args.matrix_dir / "genotypes.npy", mmap_mode="r")
    positions = np.load(args.matrix_dir / "positions.npy", mmap_mode="r")
    pairs = choose_pairs(genotypes.shape[0], args.pairs, args.pair_seed)
    blocks = [
        (start, min(start + args.block_length, args.sequence_length))
        for start in range(0, args.sequence_length, args.block_length)
    ]
    means, variances, index_map = translate_from_genotype_matrix(
        genotypes,
        positions,
        model,
        blocks=blocks,
        pivot_pairs=[tuple(map(int, pair)) for pair in pairs],
        mutation_rate=args.mutation_rate,
        device="cuda",
        batch_size=128,
        build_workers=6,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.out,
        means=means.astype(np.float32),
        variances=variances.astype(np.float32),
        index_map=index_map.astype(np.int32),
        block_starts=np.asarray([left for left, _ in blocks], dtype=np.float64),
        block_ends=np.asarray([right for _, right in blocks], dtype=np.float64),
        pairs=pairs,
        mutation_rate=np.float64(args.mutation_rate),
    )
    metadata = {
        "checkpoint": str(checkpoint),
        "checkpoint_sha256": sha256(checkpoint),
        "selected_epoch": int(float(validation["epoch"])),
        "validation_loss": float(validation["val_loss"]),
        "pairs": len(pairs),
        "blocks": len(blocks),
        "output_rows": int(len(means)),
        "output_windows": int(means.shape[1]),
        "output_sha256": sha256(args.out),
    }
    args.out.with_suffix(args.out.suffix + ".json").write_text(
        json.dumps(metadata, indent=2) + "\n"
    )
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
