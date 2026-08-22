from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from ghist_bottleneck.features import FeatureCollection
from ghist_bottleneck.infer import predict, predict_blocks
from ghist_bottleneck.simulate import SimulationConfig, generate_dataset, simulate_one
from ghist_bottleneck.train import train_ensemble


def main() -> None:
    root = Path("artifacts/scheduled_smoke")
    root.mkdir(parents=True, exist_ok=True)
    # Deliberately unrelated to the sealed public testing answer.
    smoke_generations = 11_000
    smoke_fraction = 0.31
    observed = simulate_one(
        n_haplotypes=40,
        generations=smoke_generations,
        fraction=smoke_fraction,
        ancestral_ne=15_000,
        cfg=SimulationConfig(sequence_length=250_000),
        seed=8317,
    )
    observed_path = root / "synthetic_observed.npz"
    observed.save(observed_path)
    dataset_path = root / "simulations.npz"
    generate_dataset(observed, dataset_path, n_sims=96, seed=51907, quick=True)
    checkpoint_dir = root / "checkpoints"
    train_ensemble(dataset_path, checkpoint_dir, quick=True)
    result = predict(observed_path, checkpoint_dir, root / "prediction.yaml")
    observed_blocks = [
        simulate_one(
            n_haplotypes=40,
            generations=smoke_generations,
            fraction=smoke_fraction,
            ancestral_ne=15_000,
            cfg=SimulationConfig(sequence_length=250_000),
            seed=8317 + index,
        )
        for index in range(4)
    ]
    collection = FeatureCollection(
        windows=np.asarray([block.windows for block in observed_blocks]),
        sfs=np.asarray([block.sfs for block in observed_blocks]),
        n_haplotypes=40,
        sequence_length=1_000_000,
        block_length=250_000,
        block_starts=np.arange(4) * 250_000,
        block_ends=np.arange(1, 5) * 250_000,
        n_variants=np.asarray([block.n_variants for block in observed_blocks]),
    )
    collection_path = root / "synthetic_observed_blocks.npz"
    collection.save(collection_path)
    block_result = predict_blocks(
        collection_path, checkpoint_dir, root / "prediction.blocks.yaml"
    )
    report = {"single": result, "blocks": block_result}
    (root / "result.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
