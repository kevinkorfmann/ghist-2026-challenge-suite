from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import msprime
import numpy as np
from tqdm import tqdm

from .features import FeatureBundle, FeatureCollection, LDAccumulator, summarize_variants
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
    discrete_genome: bool = False


def _recombination_model(
    cfg: SimulationConfig,
    *,
    multiplier: float,
    coefficient_of_variation: float,
    rng: np.random.Generator,
) -> float | msprime.RateMap:
    mean_rate = cfg.recombination_rate * multiplier
    if coefficient_of_variation <= 0:
        return mean_rate
    # A lognormal map is strictly positive and is renormalized so local
    # heterogeneity cannot silently change the announced chromosome mean.
    sigma = np.sqrt(np.log1p(coefficient_of_variation**2))
    rates = rng.lognormal(-0.5 * sigma**2, sigma, 32)
    rates *= mean_rate / rates.mean()
    return msprime.RateMap(
        position=np.linspace(0, cfg.sequence_length, len(rates) + 1),
        rate=rates,
    )


def _apply_spatial_suppression(
    bundle: FeatureBundle,
    *,
    intensity: float,
    rng: np.random.Generator,
) -> None:
    """Add smooth low-diversity tracts without changing the target history."""
    if intensity <= 0:
        return
    field = rng.normal(size=len(bundle.windows))
    kernel = np.array([1, 2, 3, 2, 1], dtype=float)
    kernel /= kernel.sum()
    field = np.convolve(np.pad(field, 2, mode="reflect"), kernel, mode="valid")
    field = 1 / (1 + np.exp(-field))
    retention = np.clip(1 - intensity * field, 0.25, 1.0)
    density = np.expm1(bundle.windows[:, 0])
    bundle.windows[:, 0] = np.log1p(density * retention)
    bundle.windows[:, 1] *= retention
    bundle.windows[:, 2] *= retention
    bundle.windows[:, 7] *= np.sqrt(retention)


def _draw_prior(rng: np.random.Generator, cfg: SimulationConfig) -> tuple[float, float, float]:
    generations = np.exp(rng.uniform(np.log(cfg.generation_min), np.log(cfg.generation_max)))
    fraction = rng.uniform(cfg.fraction_min, cfg.fraction_max)
    ancestral_ne = np.exp(rng.uniform(np.log(cfg.ancestral_ne_min), np.log(cfg.ancestral_ne_max)))
    return float(generations), float(fraction), float(ancestral_ne)


def _draw_nuisance(rng: np.random.Generator) -> tuple[int, float, float, float, float, float]:
    draw = rng.random()
    if draw < 0.5:
        return 0, 1.0, 1.0, 0.0, 0.0, 0.0
    scenario = 1 if draw < 0.75 else 2 if draw < 0.9 else 3
    return (
        scenario,
        float(np.exp(rng.normal(0, 0.08))),
        float(np.exp(rng.normal(0, 0.15))),
        float(rng.uniform(0.0, 0.8)) if scenario in (1, 3) else 0.0,
        float(rng.uniform(0.0, 0.06)) if scenario in (2, 3) else 0.0,
        float(rng.uniform(0.0, 0.4)) if scenario in (2, 3) else 0.0,
    )


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
    *,
    mutation_multiplier: float = 1.0,
    recombination_multiplier: float = 1.0,
    recombination_cv: float = 0.0,
    missingness: float = 0.0,
    spatial_suppression: float = 0.0,
) -> FeatureBundle:
    rng = np.random.default_rng(seed + 2)
    demography = msprime.Demography()
    demography.add_population(name="pop", initial_size=ancestral_ne * fraction)
    demography.add_population_parameters_change(
        time=generations, initial_size=ancestral_ne, population="pop"
    )
    ancestry = msprime.sim_ancestry(
        # Sample gene copies without changing the diploid population time scale.
        samples=[msprime.SampleSet(n_haplotypes, population="pop", ploidy=1)],
        ploidy=2,
        sequence_length=cfg.sequence_length,
        recombination_rate=_recombination_model(
            cfg,
            multiplier=recombination_multiplier,
            coefficient_of_variation=recombination_cv,
            rng=rng,
        ),
        demography=demography,
        random_seed=seed,
    )
    ts = msprime.sim_mutations(
        ancestry,
        rate=cfg.mutation_rate * mutation_multiplier,
        random_seed=seed + 1,
        discrete_genome=cfg.discrete_genome,
    )
    positions, ac, called, het = [], [], [], []
    ld_accumulator = LDAccumulator()
    for variant in ts.variants():
        g = variant.genotypes
        missing = rng.random(n_haplotypes) < missingness if missingness > 0 else np.zeros(
            n_haplotypes, dtype=bool
        )
        observed = ~missing
        alt = int(np.sum(g[observed] == 1))
        n_called = int(observed.sum())
        if 0 < alt < n_called:
            positions.append(variant.site.position)
            ac.append(alt)
            called.append(n_called)
            paired_n = (n_haplotypes // 2) * 2
            pairs = g[:paired_n].reshape(-1, 2)
            pair_missing = missing[:paired_n].reshape(-1, 2).any(axis=1)
            called_pairs = pairs[~pair_missing]
            het.append(
                float(np.mean(called_pairs[:, 0] != called_pairs[:, 1]))
                if len(called_pairs)
                else 0.0
            )
            dosage = pairs.sum(axis=1).astype(np.float32)
            dosage[pair_missing] = np.nan
            ld_accumulator.add(variant.site.position, dosage)
    bundle = summarize_variants(
        np.asarray(positions), np.asarray(ac), np.asarray(called), np.asarray(het),
        sequence_length=cfg.sequence_length,
        n_haplotypes=n_haplotypes,
        n_windows=cfg.n_windows,
        sfs_bins=cfg.sfs_bins,
        ld_decay=ld_accumulator.result(),
        arg_hazards=coalescence_hazard_features(ts),
    )
    _apply_spatial_suppression(bundle, intensity=spatial_suppression, rng=rng)
    return bundle


def generate_dataset(
    observed: FeatureBundle | FeatureCollection,
    out: str | Path,
    *,
    n_sims: int,
    seed: int = 20260821,
    quick: bool = False,
) -> None:
    cfg = SimulationConfig(sequence_length=250_000 if quick else 2_000_000)
    rng = np.random.default_rng(seed)
    windows, spectra, targets, nuisance = [], [], [], []
    nuisance_records = []
    for _ in tqdm(range(n_sims), desc="coalescent simulations"):
        generations, fraction, ancestral_ne = _draw_prior(rng, cfg)
        sim_seed = int(rng.integers(1, 2**31 - 2))
        (
            scenario,
            mutation_multiplier,
            recombination_multiplier,
            recombination_cv,
            missingness,
            spatial_suppression,
        ) = _draw_nuisance(rng)
        b = simulate_one(
            observed.n_haplotypes,
            generations,
            fraction,
            ancestral_ne,
            cfg,
            sim_seed,
            mutation_multiplier=mutation_multiplier,
            recombination_multiplier=recombination_multiplier,
            recombination_cv=recombination_cv,
            missingness=missingness,
            spatial_suppression=spatial_suppression,
        )
        _degrade_arg_context(b.sfs, rng)
        # Train for graceful degradation when genealogical inference is absent or
        # unreliable. The final element is the ARG-availability mask.
        if rng.random() < 0.2:
            b.sfs[-33:] = 0
        windows.append(b.windows)
        spectra.append(b.sfs)
        targets.append((np.log(generations), np.log(fraction / (1 - fraction))))
        nuisance.append(ancestral_ne)
        nuisance_records.append(
            (
                scenario,
                mutation_multiplier,
                recombination_multiplier,
                recombination_cv,
                missingness,
                spatial_suppression,
            )
        )
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        out,
        windows=np.asarray(windows, dtype=np.float32),
        sfs=np.asarray(spectra, dtype=np.float32),
        targets=np.asarray(targets, dtype=np.float32),
        ancestral_ne=np.asarray(nuisance, dtype=np.float32),
        nuisance=np.asarray(nuisance_records, dtype=np.float32),
        nuisance_names=np.asarray(
            [
                "scenario",
                "mutation_multiplier",
                "recombination_multiplier",
                "recombination_cv",
                "missingness",
                "spatial_suppression",
            ]
        ),
        n_haplotypes=np.int64(observed.n_haplotypes),
        coalescent_ploidy=np.int64(2),
        seed=np.int64(seed),
    )


def generate_chromosome_dataset(
    observed: FeatureBundle | FeatureCollection,
    out: str | Path,
    *,
    n_cells: int,
    blocks_per_cell: int = 50,
    cell_start: int = 0,
    seed: int = 20261001,
) -> None:
    """Generate 2 Mb blocks grouped into synthetic 100 Mb parameter cells."""
    cfg = SimulationConfig(sequence_length=2_000_000)
    rng = np.random.default_rng(seed)
    windows, spectra, targets, ancestral_sizes = [], [], [], []
    nuisance, cell_ids, variant_counts = [], [], []
    for local_cell in tqdm(range(n_cells), desc="chromosome parameter cells"):
        cell_id = cell_start + local_cell
        generations, fraction, ancestral_ne = _draw_prior(rng, cfg)
        nuisance_draw = _draw_nuisance(rng)
        for _ in range(blocks_per_cell):
            sim_seed = int(rng.integers(1, 2**31 - 2))
            (
                scenario,
                mutation_multiplier,
                recombination_multiplier,
                recombination_cv,
                missingness,
                spatial_suppression,
            ) = nuisance_draw
            bundle = simulate_one(
                observed.n_haplotypes,
                generations,
                fraction,
                ancestral_ne,
                cfg,
                sim_seed,
                mutation_multiplier=mutation_multiplier,
                recombination_multiplier=recombination_multiplier,
                recombination_cv=recombination_cv,
                missingness=missingness,
                spatial_suppression=spatial_suppression,
            )
            _degrade_arg_context(bundle.sfs, rng)
            if rng.random() < 0.2:
                bundle.sfs[-33:] = 0
            windows.append(bundle.windows)
            spectra.append(bundle.sfs)
            targets.append((np.log(generations), np.log(fraction / (1 - fraction))))
            ancestral_sizes.append(ancestral_ne)
            nuisance.append(nuisance_draw)
            cell_ids.append(cell_id)
            variant_counts.append(bundle.n_variants)
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        out,
        windows=np.asarray(windows, dtype=np.float32),
        sfs=np.asarray(spectra, dtype=np.float32),
        targets=np.asarray(targets, dtype=np.float32),
        ancestral_ne=np.asarray(ancestral_sizes, dtype=np.float32),
        nuisance=np.asarray(nuisance, dtype=np.float32),
        nuisance_names=np.asarray(
            [
                "scenario",
                "mutation_multiplier",
                "recombination_multiplier",
                "recombination_cv",
                "missingness",
                "spatial_suppression",
            ]
        ),
        cell_id=np.asarray(cell_ids, dtype=np.int64),
        n_variants=np.asarray(variant_counts, dtype=np.int64),
        n_haplotypes=np.int64(observed.n_haplotypes),
        coalescent_ploidy=np.int64(2),
        blocks_per_cell=np.int64(blocks_per_cell),
        seed=np.int64(seed),
    )
