from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path

import numpy as np

from ghist_bottleneck.sfs_baseline import (
    complete_block_spectra_from_vcf,
    fit_two_epoch_spectrum,
)


def _fit_replica(arguments: tuple) -> dict:
    spectrum, mutation_rate, sequence_length, folded, initial, seed = arguments
    return fit_two_epoch_spectrum(
        spectrum,
        mutation_rate=mutation_rate,
        sequence_length=sequence_length,
        restarts=2,
        seed=seed,
        folded=folded,
        initial=tuple(initial),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--fit", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--replicates", type=int, default=200)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--seed", type=int, default=83_771)
    parser.add_argument("--folded", action="store_true")
    parser.add_argument("--mutation-rate", type=float, default=3.42e-8)
    parser.add_argument("--sequence-length", type=int, default=100_000_000)
    parser.add_argument("--block-length", type=int, default=2_000_000)
    args = parser.parse_args()

    reference = json.loads(args.fit.read_text())
    if bool(reference["folded"]) != args.folded:
        raise ValueError("reference fit and bootstrap folding settings disagree")
    blocks, metadata = complete_block_spectra_from_vcf(
        args.vcf,
        sequence_length=args.sequence_length,
        block_length=args.block_length,
        polarized=bool(reference["polarized"]),
    )
    rng = np.random.default_rng(args.seed)
    tasks = []
    initial = (reference["post_decline_fraction"], reference["scaled_time_2N"])
    for replica in range(args.replicates):
        selected = rng.integers(0, len(blocks), len(blocks))
        tasks.append(
            (
                blocks[selected].sum(axis=0),
                args.mutation_rate,
                args.sequence_length,
                args.folded,
                initial,
                args.seed + replica * 7_919,
            )
        )
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        fits = list(executor.map(_fit_replica, tasks))
    values = np.asarray(
        [[fit["generations"], fit["post_decline_fraction"], fit["ancestral_effective_size"]] for fit in fits]
    )
    bounded_fraction = np.clip(values[:, 1], 1e-6, 1 - 1e-6)
    latent = np.column_stack(
        (np.log(values[:, 0]), np.log(bounded_fraction / (1 - bounded_fraction)))
    )
    result = {
        "replicates": args.replicates,
        "seed": args.seed,
        "folded": args.folded,
        "block_metadata": metadata,
        "columns": ["generations", "post_decline_fraction", "ancestral_effective_size"],
        "mean": values.mean(axis=0).tolist(),
        "median": np.median(values, axis=0).tolist(),
        "interval_95": np.quantile(values, [0.025, 0.975], axis=0).T.tolist(),
        "latent_covariance": np.cov(latent, rowvar=False).tolist(),
        "reference": {
            "generations": reference["generations"],
            "post_decline_fraction": reference["post_decline_fraction"],
            "ancestral_effective_size": reference["ancestral_effective_size"],
        },
        "values": values.tolist(),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: result[key] for key in ("replicates", "mean", "median", "interval_95")}, indent=2))


if __name__ == "__main__":
    main()
