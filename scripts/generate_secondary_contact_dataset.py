from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path

import numpy as np

from ghist_secondary_contact.demography import SecondaryContactParameters
from ghist_secondary_contact.simulate import draw_nuisance, parameter_design, simulate_block


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _simulate(arguments: tuple[int, int, int, int, float]) -> dict:
    cell_id, n_points, replicates, seed, sequence_length = arguments
    points, partitions = parameter_design(n_points, seed=seed)
    point_id, replicate = divmod(cell_id, replicates)
    ancestral, size1, size2, split, contact, migration = points[point_id]
    parameters = SecondaryContactParameters(
        ancestral_size=ancestral,
        population1_size=size1,
        population2_size=size2,
        split_time=split,
        contact_time=contact,
        migration_rate=migration,
    )
    rng = np.random.default_rng(seed + 1_000_003 * cell_id)
    nuisance = draw_nuisance(rng)
    simulation_seed = int(rng.integers(1, 2**31 - 2))
    bundle = simulate_block(
        parameters,
        nuisance,
        seed=simulation_seed,
        sequence_length=sequence_length,
        n_windows=50,
    )
    return {
        "cell_id": cell_id,
        "point_id": point_id,
        "replicate": replicate,
        "partition": int(partitions[point_id]),
        "truth": np.asarray((size1, size2, split, contact, migration)),
        "ancestral_size": ancestral,
        "nuisance": np.asarray(
            (
                nuisance.mutation_multiplier,
                nuisance.recombination_multiplier,
                nuisance.migration_asymmetry,
                *(nuisance.ancestral_epoch_factors or (np.nan,) * 3),
                *(nuisance.population1_epoch_factors or (np.nan,) * 2),
                np.nan
                if nuisance.population2_epoch_factor is None
                else nuisance.population2_epoch_factor,
                np.nan
                if nuisance.population2_contact_factor is None
                else nuisance.population2_contact_factor,
            ),
            dtype=np.float64,
        ),
        "windows": bundle.windows,
        "joint_sfs": bundle.joint_sfs,
        "global_summary": bundle.global_summary,
        "n_variants": bundle.n_variants,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--cell-start", type=int, required=True)
    parser.add_argument("--cells", type=int, default=4)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--points", type=int, default=600)
    parser.add_argument("--replicates", type=int, default=6)
    parser.add_argument("--sequence-length", type=float, default=100_000_000)
    parser.add_argument("--seed", type=int, default=20260822)
    args = parser.parse_args()
    total = args.points * args.replicates
    stop = args.cell_start + args.cells
    if args.cell_start < 0 or stop > total:
        raise ValueError(f"requested cells must lie within [0, {total})")
    work = [
        (cell_id, args.points, args.replicates, args.seed, args.sequence_length)
        for cell_id in range(args.cell_start, stop)
    ]
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        records = list(executor.map(_simulate, work))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    nuisance = np.asarray([record["nuisance"] for record in records], dtype=np.float64)
    np.savez_compressed(
        args.out,
        cell_id=np.asarray([record["cell_id"] for record in records]),
        point_id=np.asarray([record["point_id"] for record in records]),
        replicate=np.asarray([record["replicate"] for record in records]),
        partition=np.asarray([record["partition"] for record in records], dtype=np.int8),
        truth=np.asarray([record["truth"] for record in records], dtype=np.float64),
        ancestral_size=np.asarray([record["ancestral_size"] for record in records]),
        nuisance=nuisance,
        windows=np.asarray([record["windows"] for record in records], dtype=np.float32),
        joint_sfs=np.asarray([record["joint_sfs"] for record in records], dtype=np.float32),
        global_summary=np.asarray(
            [record["global_summary"] for record in records], dtype=np.float32
        ),
        n_variants=np.asarray([record["n_variants"] for record in records]),
    )
    metadata = {
        "path": str(args.out),
        "sha256": _sha256(args.out),
        "cell_start": args.cell_start,
        "cells": args.cells,
        "points": args.points,
        "replicates": args.replicates,
        "total_cells": total,
        "sequence_length": args.sequence_length,
        "seed": args.seed,
        "partition_labels": {"0": "calibration", "1": "selection", "2": "evaluation"},
    }
    args.out.with_suffix(args.out.suffix + ".json").write_text(
        json.dumps(metadata, indent=2) + "\n"
    )
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
