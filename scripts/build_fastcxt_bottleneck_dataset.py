from __future__ import annotations

import argparse
from pathlib import Path

from ghist_bottleneck.fastcxt_dataset import FastCxtSimulationConfig, build_dataset


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--cells", type=int, default=256)
    parser.add_argument("--replicas", type=int, default=2)
    parser.add_argument("--pairs", type=int, default=100)
    parser.add_argument("--workers", type=int, default=24)
    parser.add_argument("--seed", type=int, default=2_026_082_1)
    args = parser.parse_args()
    config = FastCxtSimulationConfig(num_pairs=args.pairs, master_seed=args.seed)
    build_dataset(
        args.out,
        n_cells=args.cells,
        replicas=args.replicas,
        workers=args.workers,
        config=config,
    )


if __name__ == "__main__":
    main()
