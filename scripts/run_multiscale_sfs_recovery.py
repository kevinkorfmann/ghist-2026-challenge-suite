from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path

import msprime
import numpy as np
from scipy.stats import qmc

from ghist_bottleneck.sfs_baseline import fit_two_epoch_spectrum
from scripts.bottleneck_consumer_guard import (
    ROOT, bundle_metadata, capture_inputs, require_new_outputs, require_unchanged,
    source_manifest,
)


BROAD_POINTS = 100
FOCUSED_POINTS = 50
REPLICATES_PER_POINT = 6
TOTAL_POINTS = BROAD_POINTS + FOCUSED_POINTS
TOTAL_CELLS = TOTAL_POINTS * REPLICATES_PER_POINT


def _scale_design(unit: np.ndarray, bounds: tuple[tuple[float, float], ...]) -> np.ndarray:
    result = np.empty_like(unit)
    for coordinate, (lower, upper) in enumerate(bounds):
        if coordinate in (0, 2):
            result[:, coordinate] = np.exp(
                np.log(lower) + unit[:, coordinate] * (np.log(upper) - np.log(lower))
            )
        else:
            result[:, coordinate] = lower + unit[:, coordinate] * (upper - lower)
    return result


def parameter_design(seed: int = 20261117) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    broad = _scale_design(
        qmc.LatinHypercube(3, seed=seed).random(BROAD_POINTS),
        ((200.0, 30000.0), (0.08, 0.95), (5000.0, 50000.0)),
    )
    focused = _scale_design(
        qmc.LatinHypercube(3, seed=seed + 1).random(FOCUSED_POINTS),
        ((3000.0, 6500.0), (0.42, 0.62), (15000.0, 30000.0)),
    )
    points = np.vstack((broad, focused))
    domain = np.concatenate(
        (np.zeros(BROAD_POINTS, dtype=np.int8), np.ones(FOCUSED_POINTS, dtype=np.int8))
    )
    split = np.empty(TOTAL_POINTS, dtype=np.int8)
    rng = np.random.default_rng(seed + 2)
    for ids, calibration, selection in (
        (np.arange(BROAD_POINTS), 67, 17),
        (np.arange(BROAD_POINTS, TOTAL_POINTS), 33, 8),
    ):
        shuffled = rng.permutation(ids)
        split[shuffled[:calibration]] = 0
        split[shuffled[calibration : calibration + selection]] = 1
        split[shuffled[calibration + selection :]] = 2
    return points, domain, split


PARAMETER_POINTS, POINT_DOMAIN, POINT_SPLIT = parameter_design()


def simulation_design(*, seed, blocks, block_length, mutation_rate, restarts):
    if any(type(value) is not int or value < 1 for value in (seed, blocks, restarts)):
        raise ValueError("Simulation seed, blocks and restarts must be positive integers")
    if not np.isfinite([block_length, mutation_rate]).all() or min(block_length, mutation_rate) <= 0:
        raise ValueError("Block length and mutation rate must be finite and positive")
    return {
        "schema": "multiscale_independent_blocks_v1", "seed": seed,
        "blocks": blocks, "block_length": float(block_length),
        "sequence_length": blocks * float(block_length), "mutation_rate": float(mutation_rate),
        "recombination_rate": 1.87e-8, "restarts": restarts,
        "sample_haplotypes": 114, "coalescent_ploidy": 2,
        "parameter_design_seed": 20261117,
        "parameter_design_sha256": hashlib.sha256(PARAMETER_POINTS.astype("<f8").tobytes()
            + POINT_DOMAIN.tobytes() + POINT_SPLIT.tobytes()).hexdigest(),
        "simulation_seed_rule": "default_rng(seed + cell_id * 104729), consecutive ancestry/mutation draws per block",
        "optimizer_seed_rule": "seed + 900000 + cell_id",
        "linkage": "independent ancestry simulation per block; not one linked chromosome",
    }


def require_simulation_design(array, fits, context):
    design = context.get("simulation_design")
    if not isinstance(design, dict):
        raise ValueError("Recovery shard lacks its simulation/block design")
    keys = ("seed", "blocks", "block_length", "mutation_rate", "restarts")
    if any(key not in design for key in keys) or design != simulation_design(**{key: design[key] for key in keys}):
        raise ValueError("Recovery shard simulation/block design is invalid or stale")
    expected_data = {"folded": True, "mutation_rate": design["mutation_rate"],
                     "sequence_length": design["sequence_length"], "sample_haplotypes": 114}
    if context.get("fit_data_contract") != expected_data:
        raise ValueError("Recovery fit exposure or sampling differs from its simulation design")
    for cell, fit in zip(array["cell_id"], fits, strict=True):
        if fit.get("optimizer_seed") != design["seed"] + 900000 + int(cell) or fit.get("restarts") != design["restarts"]:
            raise ValueError("Recovery fit seed/restarts differ from the declared design")
    return design


def multiscale_parameters(cell_id: int) -> tuple[int, int, int, int, float, float, float]:
    if not 0 <= cell_id < TOTAL_CELLS:
        raise ValueError(f"cell_id must be in [0, {TOTAL_CELLS})")
    point_id, replicate = divmod(cell_id, REPLICATES_PER_POINT)
    generations, fraction, ancestral_size = PARAMETER_POINTS[point_id]
    return (
        int(point_id),
        int(replicate),
        int(POINT_DOMAIN[point_id]),
        int(POINT_SPLIT[point_id]),
        float(generations),
        float(fraction),
        float(ancestral_size),
    )


def simulate_cell(arguments: tuple[int, int, int, float, float, int]) -> dict:
    cell_id, seed, blocks, block_length, mutation_rate, restarts = arguments
    point_id, replicate, domain, split, generations, fraction, ancestral_size = (
        multiscale_parameters(cell_id)
    )
    rng = np.random.default_rng(seed + cell_id * 104729)
    demography = msprime.Demography()
    demography.add_population(name="pop", initial_size=ancestral_size * fraction)
    demography.add_population_parameters_change(
        time=generations, initial_size=ancestral_size, population="pop"
    )
    spectrum = np.zeros(115, dtype=np.float64)
    for _ in range(blocks):
        ancestry_seed = int(rng.integers(1, 2**31 - 2))
        mutation_seed = int(rng.integers(1, 2**31 - 2))
        ancestry = msprime.sim_ancestry(
            samples=[msprime.SampleSet(114, population="pop", ploidy=1)],
            ploidy=2,
            sequence_length=block_length,
            recombination_rate=1.87e-8,
            demography=demography,
            random_seed=ancestry_seed,
        )
        mutated = msprime.sim_mutations(
            ancestry,
            rate=mutation_rate,
            random_seed=mutation_seed,
            model=msprime.BinaryMutationModel(),
            discrete_genome=False,
        )
        spectrum += mutated.allele_frequency_spectrum(
            mode="site", polarised=True, span_normalise=False
        )
    fit = fit_two_epoch_spectrum(
        spectrum,
        mutation_rate=mutation_rate,
        sequence_length=blocks * block_length,
        restarts=restarts,
        seed=seed + 900_000 + cell_id,
        folded=True,
    )
    return {
        "cell_id": cell_id,
        "point_id": point_id,
        "replicate": replicate,
        "domain": domain,
        "split": split,
        "truth": [generations, fraction, ancestral_size],
        "estimate": [
            fit["generations"],
            fit["post_decline_fraction"],
            fit["ancestral_effective_size"],
        ],
        "scaled_time_2N": fit["scaled_time_2N"],
        "theta": fit["theta"],
        "multinomial_log_likelihood": fit["multinomial_log_likelihood"],
        "segregating_sites": int(spectrum[1:-1].sum()),
        "fit": fit,
    }


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def recorded_cell(arguments):
    try:
        return simulate_cell(arguments)
    except Exception as error:
        return {"cell_id": arguments[0], "failure": f"{type(error).__name__}: {error}"}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--cell-start", type=int, required=True)
    parser.add_argument("--cells", type=int, default=10)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--blocks", type=int, default=50)
    parser.add_argument("--block-length", type=float, default=2_000_000)
    parser.add_argument("--mutation-rate", type=float, default=3.42e-8)
    parser.add_argument("--restarts", type=int, default=8)
    parser.add_argument("--seed", type=int, default=20261117)
    args = parser.parse_args()
    design = simulation_design(seed=args.seed, blocks=args.blocks, block_length=args.block_length,
                               mutation_rate=args.mutation_rate, restarts=args.restarts)
    diagnostics = args.out.with_suffix(args.out.suffix + ".fits.json")
    require_new_outputs(args.out, args.out.with_suffix(args.out.suffix + ".json"), diagnostics)
    source_hashes = capture_inputs(*(ROOT / name for name in source_manifest(
        "scripts/run_multiscale_sfs_recovery.py")))
    stop = args.cell_start + args.cells
    if args.cell_start < 0 or stop > TOTAL_CELLS:
        raise ValueError(f"requested cells must lie within [0, {TOTAL_CELLS})")
    work = [
        (cell_id, args.seed, args.blocks, args.block_length, args.mutation_rate, args.restarts)
        for cell_id in range(args.cell_start, stop)
    ]
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        records = list(executor.map(recorded_cell, work))
    # Keep numerical failures inspectable; no incomplete/unsupported shard can
    # enter calibration merely because its physical estimates are finite.
    diagnostics.parent.mkdir(parents=True, exist_ok=True)
    diagnostics.write_text(json.dumps({"planned_cells": [task[0] for task in work],
                                      "cell_fits": records}, indent=2, allow_nan=False) + "\n")
    if any("failure" in record for record in records):
        raise ValueError(f"Recovery cells failed; diagnostics retained in {diagnostics}")
    provenance = bundle_metadata([record["fit"] for record in records],
                                 producer="scripts/run_multiscale_sfs_recovery.py", simulation_design=design)
    require_unchanged(source_hashes)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.out,
        cell_id=np.asarray([record["cell_id"] for record in records], dtype=np.int64),
        point_id=np.asarray([record["point_id"] for record in records], dtype=np.int64),
        replicate=np.asarray([record["replicate"] for record in records], dtype=np.int64),
        domain=np.asarray([record["domain"] for record in records], dtype=np.int8),
        split=np.asarray([record["split"] for record in records], dtype=np.int8),
        truth=np.asarray([record["truth"] for record in records], dtype=np.float64),
        estimate=np.asarray([record["estimate"] for record in records], dtype=np.float64),
        scaled_time_2N=np.asarray([record["scaled_time_2N"] for record in records]),
        theta=np.asarray([record["theta"] for record in records]),
        multinomial_log_likelihood=np.asarray(
            [record["multinomial_log_likelihood"] for record in records]
        ),
        segregating_sites=np.asarray([record["segregating_sites"] for record in records]),
        coalescent_ploidy=np.int64(2),
        **provenance,
    )
    metadata = {
        "output": str(args.out),
        "sha256": sha256(args.out),
        "fit_diagnostics_sha256": sha256(diagnostics),
        "optimizer_contract": json.loads(str(provenance["optimizer_contract_json"])),
        "source_sha256": source_manifest("scripts/run_multiscale_sfs_recovery.py"),
        "simulation_design": design,
        "scientifically_validated": False,
        "release_ready": False,
        "cell_start": args.cell_start,
        "cells": args.cells,
        "workers": args.workers,
        "blocks": args.blocks,
        "block_length": args.block_length,
        "sequence_length": args.blocks * args.block_length,
        "mutation_rate": args.mutation_rate,
        "coalescent_ploidy": 2,
        "restarts": args.restarts,
        "seed": args.seed,
        "broad_points": BROAD_POINTS,
        "focused_points": FOCUSED_POINTS,
        "replicates_per_point": REPLICATES_PER_POINT,
        "total_cells": TOTAL_CELLS,
        "domain_labels": {"0": "broad", "1": "focused"},
        "split_labels": {"0": "calibration", "1": "selection", "2": "evaluation"},
    }
    args.out.with_suffix(args.out.suffix + ".json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
