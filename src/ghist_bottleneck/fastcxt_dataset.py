from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import msprime
import numpy as np
import tskit


@dataclass(frozen=True)
class FastCxtSimulationConfig:
    n_haplotypes: int = 114
    sequence_length: int = 1_000_000
    window_size: int = 2_000
    mutation_rate: float = 3.42e-8
    recombination_rate: float = 1.87e-8
    generations_min: float = 200.0
    generations_max: float = 30_000.0
    fraction_min: float = 0.08
    fraction_max: float = 0.95
    ancestral_ne_min: float = 5_000.0
    ancestral_ne_max: float = 50_000.0
    num_pairs: int = 100
    validation_fraction: float = 0.2
    master_seed: int = 2_026_082_1


def grouped_splits(n_cells: int, validation_fraction: float, seed: int) -> np.ndarray:
    if n_cells < 2 or not 0 < validation_fraction < 1:
        raise ValueError("at least two cells and a validation fraction in (0, 1) are required")
    n_validation = max(1, min(n_cells - 1, round(n_cells * validation_fraction)))
    order = np.random.default_rng(seed).permutation(n_cells)
    splits = np.full(n_cells, "train", dtype="U5")
    splits[order[:n_validation]] = "test"
    return splits


def choose_pairs(n_samples: int, num_pairs: int, seed: int) -> np.ndarray:
    pairs = np.asarray(
        [(left, right) for left in range(n_samples) for right in range(left + 1, n_samples)],
        dtype=np.int32,
    )
    if len(pairs) <= num_pairs:
        return pairs
    selected = np.random.default_rng(seed).choice(len(pairs), num_pairs, replace=False)
    return pairs[selected]


def windowed_tmrca(
    ts: tskit.TreeSequence,
    sample_a: int,
    sample_b: int,
    window_size: int,
) -> np.ndarray:
    n_windows = int(np.ceil(ts.sequence_length / window_size))
    total = np.zeros(n_windows, dtype=np.float64)
    for tree in ts.trees():
        node = tree.mrca(sample_a, sample_b)
        if node == tskit.NULL:
            continue
        left, right = tree.interval
        first = int(left // window_size)
        last = min(n_windows - 1, int(np.nextafter(right, left) // window_size))
        for window in range(first, last + 1):
            window_left = window * window_size
            window_right = min((window + 1) * window_size, ts.sequence_length)
            overlap = min(right, window_right) - max(left, window_left)
            if overlap > 0:
                total[window] += overlap * tree.time(node)
    widths = np.minimum(
        window_size,
        ts.sequence_length - np.arange(n_windows, dtype=np.float64) * window_size,
    )
    return total / widths


def _log_uniform(rng: np.random.Generator, low: float, high: float) -> float:
    return float(np.exp(rng.uniform(np.log(low), np.log(high))))


def _simulate_cell_replica(args: tuple) -> dict:
    cell, replica, split, out_root, cfg = args
    cfg = FastCxtSimulationConfig(**cfg)
    cell_rng = np.random.default_rng(cfg.master_seed + cell * 1_000_003)
    generations = _log_uniform(cell_rng, cfg.generations_min, cfg.generations_max)
    fraction = float(cell_rng.uniform(cfg.fraction_min, cfg.fraction_max))
    ancestral_ne = _log_uniform(cell_rng, cfg.ancestral_ne_min, cfg.ancestral_ne_max)

    rng = np.random.default_rng(cfg.master_seed + cell * 1_000_003 + replica * 10_007 + 17)
    scenario = "nominal" if replica == 0 else "rate_and_coverage"
    mutation_multiplier = 1.0 if replica == 0 else float(np.exp(rng.normal(0, 0.08)))
    recombination_multiplier = 1.0 if replica == 0 else float(np.exp(rng.normal(0, 0.15)))
    recombination_cv = 0.0 if replica == 0 else float(rng.uniform(0.0, 0.8))
    site_dropout = 0.0 if replica == 0 else float(rng.uniform(0.0, 0.08))
    tract_suppression = 0.0 if replica == 0 else float(rng.uniform(0.0, 0.4))

    demography = msprime.Demography()
    demography.add_population(name="pop", initial_size=ancestral_ne * fraction)
    demography.add_population_parameters_change(
        time=generations, initial_size=ancestral_ne, population="pop"
    )
    recombination: float | msprime.RateMap
    if recombination_cv > 0:
        sigma = np.sqrt(np.log1p(recombination_cv**2))
        rates = rng.lognormal(-0.5 * sigma**2, sigma, 32)
        rates *= cfg.recombination_rate * recombination_multiplier / rates.mean()
        recombination = msprime.RateMap(
            position=np.linspace(0, cfg.sequence_length, len(rates) + 1), rate=rates
        )
    else:
        recombination = cfg.recombination_rate * recombination_multiplier
    ancestry = msprime.sim_ancestry(
        samples={"pop": cfg.n_haplotypes},
        ploidy=1,
        demography=demography,
        sequence_length=cfg.sequence_length,
        recombination_rate=recombination,
        random_seed=int(rng.integers(1, 2**31 - 1)),
    )
    mutated = msprime.sim_mutations(
        ancestry,
        rate=cfg.mutation_rate * mutation_multiplier,
        discrete_genome=False,
        random_seed=int(rng.integers(1, 2**31 - 1)),
    )
    genotypes = mutated.genotype_matrix().T.astype(np.int8, copy=False)
    positions = mutated.tables.sites.position.copy()
    frequency = genotypes.sum(axis=0)
    keep = (frequency > 0) & (frequency < cfg.n_haplotypes)
    if site_dropout > 0:
        keep &= rng.random(len(keep)) >= site_dropout
    if tract_suppression > 0:
        center = rng.uniform(0, cfg.sequence_length)
        scale = rng.uniform(0.08, 0.25) * cfg.sequence_length
        field = np.exp(-0.5 * np.square((positions - center) / scale))
        keep &= rng.random(len(keep)) >= tract_suppression * field
    genotypes, positions = genotypes[:, keep], positions[keep]

    pairs = choose_pairs(cfg.n_haplotypes, cfg.num_pairs, int(rng.integers(1, 2**31 - 1)))
    targets = np.empty((len(pairs), int(np.ceil(cfg.sequence_length / cfg.window_size))), np.float16)
    for index, (left, right) in enumerate(pairs):
        value = windowed_tmrca(ancestry, int(left), int(right), cfg.window_size)
        targets[index] = np.log(np.clip(value, 1e-10, None)).astype(np.float16)

    directory = Path(out_root) / split / f"cell_{cell:05d}" / f"replica_{replica:02d}"
    directory.mkdir(parents=True, exist_ok=False)
    np.save(directory / "genotypes.npy", genotypes)
    np.save(directory / "positions.npy", positions)
    np.save(directory / "pairs.npy", pairs)
    np.save(directory / "y.npy", targets)
    metadata = {
        "storage_mode": "lazy_sfs",
        "cell": cell,
        "replica": replica,
        "split": split,
        "scenario": scenario,
        "num_samples": cfg.n_haplotypes,
        "sequence_length": cfg.sequence_length,
        "window_size": cfg.window_size,
        "mutation_rate": cfg.mutation_rate * mutation_multiplier,
        "mutation_multiplier": mutation_multiplier,
        "recombination_multiplier": recombination_multiplier,
        "recombination_cv": recombination_cv,
        "site_dropout": site_dropout,
        "tract_suppression": tract_suppression,
        "generations": generations,
        "fraction": fraction,
        "ancestral_ne": ancestral_ne,
        "num_sites": int(genotypes.shape[1]),
        "num_pairs": int(len(pairs)),
    }
    (directory / "meta.json").write_text(json.dumps(metadata, indent=2) + "\n")
    return metadata


def build_dataset(
    out_root: str | Path,
    *,
    n_cells: int,
    replicas: int,
    workers: int,
    config: FastCxtSimulationConfig,
) -> None:
    out_root = Path(out_root)
    out_root.mkdir(parents=True, exist_ok=False)
    splits = grouped_splits(n_cells, config.validation_fraction, config.master_seed)
    arguments = [
        (cell, replica, str(splits[cell]), str(out_root), asdict(config))
        for cell in range(n_cells)
        for replica in range(replicas)
    ]
    with ProcessPoolExecutor(max_workers=workers) as executor:
        records = list(executor.map(_simulate_cell_replica, arguments))
    manifest = {
        "config": asdict(config),
        "n_cells": n_cells,
        "replicas": replicas,
        "grouped_split": True,
        "train_cells": int(np.sum(splits == "train")),
        "test_cells": int(np.sum(splits == "test")),
        "records": records,
    }
    (out_root / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
