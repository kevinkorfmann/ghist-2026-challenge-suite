from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path

import numpy as np

from ghist_secondary_contact.simulate import parameter_design
from ghist_secondary_contact.spectrum import secondary_contact_spectrum


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _one_point(arguments: tuple[int, np.ndarray, int, float, float, int]) -> dict:
    point_id, point, replicates, mutation_rate, sequence_length, seed = arguments
    ancestral, size1, size2, split, contact, migration = point
    scaled = (
        size1 / ancestral,
        size2 / ancestral,
        (split - contact) / (2 * ancestral),
        contact / (2 * ancestral),
        2 * ancestral * migration,
    )
    expected = np.asarray(secondary_contact_spectrum(scaled, (86, 48)), dtype=np.float64)
    expected = np.clip(expected, 0.0, None)
    rng = np.random.default_rng(seed + point_id * 1_000_003)
    spectra, nuisance, counts = [], [], []
    for replicate in range(replicates):
        if replicate == 0 or rng.random() < 0.5:
            mutation_multiplier = 1.0
            orientation_error = 0.0
        else:
            mutation_multiplier = float(
                np.exp(rng.uniform(np.log(0.85), np.log(1.15)))
            )
            orientation_error = float(rng.uniform(0.0, 0.03))
        oriented = (1 - orientation_error) * expected + orientation_error * expected[::-1, ::-1]
        theta = 4 * ancestral * mutation_rate * mutation_multiplier * sequence_length
        observed = rng.poisson(theta * oriented).astype(np.float64)
        counts.append(int(observed.sum()))
        transformed = np.log1p(observed)
        transformed /= np.linalg.norm(transformed) + 1e-12
        spectra.append(transformed)
        nuisance.append((mutation_multiplier, orientation_error))
    return {
        "point_id": point_id,
        "truth": np.asarray((size1, size2, split, contact, migration)),
        "ancestral_size": ancestral,
        "spectra": np.asarray(spectra),
        "nuisance": np.asarray(nuisance),
        "variant_counts": np.asarray(counts),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--point-start", type=int, required=True)
    parser.add_argument("--point-count", type=int, default=20)
    parser.add_argument("--points", type=int, default=2000)
    parser.add_argument("--replicates", type=int, default=6)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--mutation-rate", type=float, default=8.4e-8)
    parser.add_argument("--sequence-length", type=float, default=100_000_000)
    parser.add_argument("--seed", type=int, default=20260822)
    args = parser.parse_args()
    stop = args.point_start + args.point_count
    if args.point_start < 0 or stop > args.points:
        raise ValueError(f"requested points must lie within [0, {args.points})")
    points, partitions = parameter_design(args.points, seed=args.seed)
    work = [
        (
            point_id,
            points[point_id],
            args.replicates,
            args.mutation_rate,
            args.sequence_length,
            args.seed,
        )
        for point_id in range(args.point_start, stop)
    ]
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        records = list(executor.map(_one_point, work))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.out,
        point_id=np.repeat([record["point_id"] for record in records], args.replicates),
        replicate=np.tile(np.arange(args.replicates), len(records)),
        partition=np.repeat(
            [partitions[record["point_id"]] for record in records], args.replicates
        ).astype(np.int8),
        truth=np.repeat(
            np.asarray([record["truth"] for record in records]), args.replicates, axis=0
        ),
        ancestral_size=np.repeat(
            [record["ancestral_size"] for record in records], args.replicates
        ),
        joint_sfs=np.concatenate([record["spectra"] for record in records]).astype(np.float32),
        nuisance=np.concatenate([record["nuisance"] for record in records]).astype(np.float32),
        n_variants=np.concatenate([record["variant_counts"] for record in records]),
    )
    metadata = {
        "path": str(args.out),
        "sha256": _sha256(args.out),
        "point_start": args.point_start,
        "point_count": args.point_count,
        "points": args.points,
        "replicates": args.replicates,
        "sequence_length": args.sequence_length,
        "mutation_rate": args.mutation_rate,
        "seed": args.seed,
    }
    args.out.with_suffix(args.out.suffix + ".json").write_text(
        json.dumps(metadata, indent=2) + "\n"
    )
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
