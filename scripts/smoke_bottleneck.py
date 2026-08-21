from __future__ import annotations

import json
from pathlib import Path

from ghist_bottleneck.infer import predict
from ghist_bottleneck.simulate import SimulationConfig, generate_dataset, simulate_one
from ghist_bottleneck.train import train_ensemble


def main() -> None:
    root = Path("artifacts/scheduled_smoke")
    root.mkdir(parents=True, exist_ok=True)
    observed = simulate_one(
        n_haplotypes=40,
        generations=4757,
        fraction=0.52,
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
    (root / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

