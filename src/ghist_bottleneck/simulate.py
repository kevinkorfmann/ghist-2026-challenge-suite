from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import msprime
import numpy as np
from tqdm import tqdm

from .features import FeatureBundle, LDAccumulator, summarize_variants
from .arg_features import coalescence_hazard_features


@dataclass(frozen=True)
class SimulationConfig:
    mutation_rate: float = 3.42e-8
    recombination_rate: float = 1.87e-8
    sequence_length: float = 2_000_000
    n_windows: int = 64
    sfs_bins: int = 64
    generation_min: float = 200
    generation_max: float = 30_000
    fraction_min: float = 0.08
    fraction_max: float = 0.95
    ancestral_ne_min: float = 5_000
    ancestral_ne_max: float = 50_000


def _draw_prior(rng: np.random.Generator, cfg: SimulationConfig) -> tuple[float, float, float]:
    generations = np.exp(rng.uniform(np.log(cfg.generation_min), np.log(cfg.generation_max)))
    fraction = rng.uniform(cfg.fraction_min, cfg.fraction_max)
    ancestral_ne = np.exp(rng.uniform(np.log(cfg.ancestral_ne_min), np.log(cfg.ancestral_ne_max)))
    return float(generations), float(fraction), float(ancestral_ne)


def _degrade_arg_context(context: np.ndarray, rng: np.random.Generator) -> None:
    """Mimic smoothing, time error, and uncertainty from genealogy inference."""
    hazards = context[-33:-1]
    if rng.random() < 0.7:
        sigma = rng.uniform(0.5, 2.0)
        radius = max(1, int(np.ceil(3 * sigma)))
        x = np.arange(-radius, radius + 1)
        kernel = np.exp(-0.5 * np.square(x / sigma))
        kernel /= kernel.sum()
        padded = np.pad(hazards, radius, mode="edge")
        hazards[:] = np.convolve(padded, kernel, mode="valid")
        shift = rng.uniform(-1.5, 1.5)
        grid = np.arange(len(hazards), dtype=float)
        hazards[:] = np.interp(grid + shift, grid, hazards, left=hazards[0], right=hazards[-1])
        logit = np.log(np.clip(hazards, 1e-5, 1 - 1e-5) / np.clip(1 - hazards, 1e-5, 1))
        hazards[:] = 1 / (1 + np.exp(-(logit + rng.normal(0, 0.15, len(logit)))))


def simulate_one(
    n_haplotypes: int,
    generations: float,
    fraction: float,
    ancestral_ne: float,
    cfg: SimulationConfig,
    seed: int,
) -> FeatureBundle:
    demography = msprime.Demography()
    demography.add_population(name="pop", initial_size=ancestral_ne * fraction)
    demography.add_population_parameters_change(
        time=generations, initial_size=ancestral_ne, population="pop"
    )
    ancestry = msprime.sim_ancestry(
        samples={"pop": n_haplotypes},
        ploidy=1,
        sequence_length=cfg.sequence_length,
        recombination_rate=cfg.recombination_rate,
        demography=demography,
        random_seed=seed,
    )
    ts = msprime.sim_mutations(ancestry, rate=cfg.mutation_rate, random_seed=seed + 1)
    positions, ac, called, het = [], [], [], []
    ld_accumulator = LDAccumulator()
    for variant in ts.variants():
        g = variant.genotypes
        alt = int(np.sum(g == 1))
        if 0 < alt < n_haplotypes:
            positions.append(variant.site.position)
            ac.append(alt)
            called.append(n_haplotypes)
            # Random pairing is sufficient for a polarization-free heterozygosity summary.
            pairs = g[: (n_haplotypes // 2) * 2].reshape(-1, 2)
            het.append(float(np.mean(pairs[:, 0] != pairs[:, 1])))
            dosage = pairs.sum(axis=1).astype(np.float32)
            ld_accumulator.add(variant.site.position, dosage)
    return summarize_variants(
        np.asarray(positions), np.asarray(ac), np.asarray(called), np.asarray(het),
        sequence_length=cfg.sequence_length,
        n_haplotypes=n_haplotypes,
        n_windows=cfg.n_windows,
        sfs_bins=cfg.sfs_bins,
        ld_decay=ld_accumulator.result(),
        arg_hazards=coalescence_hazard_features(ts),
    )


def generate_dataset(
    observed: FeatureBundle,
    out: str | Path,
    *,
    n_sims: int,
    seed: int = 20260821,
    quick: bool = False,
) -> None:
    cfg = SimulationConfig(sequence_length=250_000 if quick else 2_000_000)
    rng = np.random.default_rng(seed)
    windows, spectra, targets, nuisance = [], [], [], []
    for _ in tqdm(range(n_sims), desc="coalescent simulations"):
        generations, fraction, ancestral_ne = _draw_prior(rng, cfg)
        sim_seed = int(rng.integers(1, 2**31 - 2))
        b = simulate_one(observed.n_haplotypes, generations, fraction, ancestral_ne, cfg, sim_seed)
        _degrade_arg_context(b.sfs, rng)
        # Train for graceful degradation when genealogical inference is absent or
        # unreliable. The final element is the ARG-availability mask.
        if rng.random() < 0.2:
            b.sfs[-33:] = 0
        windows.append(b.windows)
        spectra.append(b.sfs)
        targets.append((np.log(generations), np.log(fraction / (1 - fraction))))
        nuisance.append(ancestral_ne)
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        out,
        windows=np.asarray(windows, dtype=np.float32),
        sfs=np.asarray(spectra, dtype=np.float32),
        targets=np.asarray(targets, dtype=np.float32),
        ancestral_ne=np.asarray(nuisance, dtype=np.float32),
        n_haplotypes=np.int64(observed.n_haplotypes),
        seed=np.int64(seed),
    )
