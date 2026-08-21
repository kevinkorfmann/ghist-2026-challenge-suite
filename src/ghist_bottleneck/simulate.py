from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import msprime
import numpy as np
from tqdm import tqdm

from .features import FeatureBundle, summarize_variants
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
    positions, ac, called, het, ld_genotypes = [], [], [], [], []
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
            ld_genotypes.append(pairs.sum(axis=1).astype(np.float32))
    return summarize_variants(
        np.asarray(positions), np.asarray(ac), np.asarray(called), np.asarray(het),
        sequence_length=cfg.sequence_length,
        n_haplotypes=n_haplotypes,
        n_windows=cfg.n_windows,
        sfs_bins=cfg.sfs_bins,
        ld_positions=np.asarray(positions),
        ld_genotypes=np.asarray(ld_genotypes),
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
