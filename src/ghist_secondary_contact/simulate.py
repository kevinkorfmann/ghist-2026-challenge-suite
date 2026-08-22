from __future__ import annotations

from dataclasses import dataclass

import msprime
import numpy as np
from scipy.stats import qmc

from .demography import SecondaryContactParameters, SizeChange, build_demography
from .features import SecondaryFeatureBundle, summarize_tree_sequence


@dataclass(frozen=True)
class SimulationNuisance:
    mutation_multiplier: float = 1.0
    recombination_multiplier: float = 1.0
    migration_asymmetry: float = 1.0
    ancestral_epoch_factors: tuple[float, float, float] = ()
    population1_epoch_factors: tuple[float, float] = ()
    population2_epoch_factor: float | None = None
    population2_contact_factor: float | None = None


def _grouped_partition(n_points: int, seed: int) -> np.ndarray:
    partition = np.full(n_points, 2, dtype=np.int8)
    shuffled = np.random.default_rng(seed + 1).permutation(n_points)
    calibration_end = round(0.70 * n_points)
    selection_end = calibration_end + round(0.15 * n_points)
    partition[shuffled[:calibration_end]] = 0
    partition[shuffled[calibration_end:selection_end]] = 1
    return partition


def parameter_design(n_points: int, *, seed: int = 20260822) -> tuple[np.ndarray, np.ndarray]:
    """Create point-grouped, space-filling targets and validation partitions."""
    if n_points < 20:
        raise ValueError("at least twenty points are required for grouped validation")
    unit = qmc.LatinHypercube(6, seed=seed).random(n_points)
    log_bounds = (
        (1_000.0, 1_000_000.0),
        (1_000.0, 1_000_000.0),
        (1_000.0, 1_000_000.0),
        (100.0, 300_000.0),
        (1e-8, 0.05),
    )
    values = np.empty((n_points, 6), dtype=np.float64)
    for output, coordinate in enumerate((0, 1, 2, 3, 5)):
        lower, upper = log_bounds[output]
        values[:, coordinate] = np.exp(
            np.log(lower) + unit[:, coordinate] * np.log(upper / lower)
        )
    lower_logit = np.log(0.002 / 0.998)
    upper_logit = np.log(0.98 / 0.02)
    fraction_logit = lower_logit + unit[:, 4] * (upper_logit - lower_logit)
    fraction = 1.0 / (1.0 + np.exp(-fraction_logit))
    values[:, 4] = values[:, 3] * fraction
    return values, _grouped_partition(n_points, seed)


def focused_parameter_design(
    n_points: int,
    center: np.ndarray,
    *,
    seed: int = 20260822,
    include_full_prior: bool = True,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Create core, guard-shell, and broad points around a mechanistic fit.

    Columns follow ancestral size, mainland size, island size, split time,
    contact time, and migration rate. Domain labels 0, 1, and 2 denote core,
    guard, and broad points respectively.
    """
    if n_points < 40:
        raise ValueError("at least forty points are required for focused validation")
    center = np.asarray(center, dtype=np.float64)
    if center.shape != (6,) or np.any(center <= 0) or center[4] >= center[3]:
        raise ValueError("invalid six-parameter focused-design center")
    unit = qmc.LatinHypercube(6, seed=seed).random(n_points)
    domain = np.full(n_points, 2, dtype=np.int8)
    order = np.random.default_rng(seed + 2).permutation(n_points)
    core_end = round(0.70 * n_points)
    guard_end = core_end + round(0.20 * n_points)
    domain[order[:core_end]] = 0
    domain[order[core_end:guard_end]] = 1
    values = np.empty((n_points, 6), dtype=np.float64)
    log_center = np.log(center[[0, 1, 2, 3, 5]])
    center_fraction = center[4] / center[3]
    center_logit = np.log(center_fraction) - np.log1p(-center_fraction)
    # Dense precision region, followed by a deliberately much wider shell.
    half_widths = {
        0: np.log(np.asarray((2.5, 3.0, 3.0, 3.0, 15.0))),
        1: np.log(np.asarray((8.0, 10.0, 10.0, 8.0, 100.0))),
        2: np.log(np.asarray((10.0, 10.0, 4.0, 10.0, 1_000.0))),
    }
    fraction_widths = {0: 2.0, 1: 4.0, 2: 6.0}
    coordinates = (0, 1, 2, 3, 5)
    broad_values, _ = parameter_design(n_points, seed=seed + 10_000)
    for row in range(n_points):
        if domain[row] == 2 and include_full_prior:
            values[row] = broad_values[row]
            continue
        width = half_widths[int(domain[row])]
        logged = log_center + (2.0 * unit[row, coordinates] - 1.0) * width
        values[row, coordinates] = np.exp(logged)
        contact_logit = center_logit + (2.0 * unit[row, 4] - 1.0) * fraction_widths[
            int(domain[row])
        ]
        fraction = 1.0 / (1.0 + np.exp(-contact_logit))
        values[row, 4] = values[row, 3] * fraction
    if not include_full_prior:
        values[:, :3] = np.clip(values[:, :3], 1_000.0, 250_000.0)
        values[:, 3] = np.clip(values[:, 3], 100.0, 100_000.0)
        values[:, 4] = np.minimum(values[:, 4], values[:, 3] * (1.0 - 1e-6))
        values[:, 5] = np.clip(values[:, 5], 1e-8, 0.05)
    return values, _grouped_partition(n_points, seed), domain


def draw_nuisance(rng: np.random.Generator) -> SimulationNuisance:
    """Sample nominal or generator-misspecified simulation conditions."""
    if rng.random() < 0.5:
        return SimulationNuisance()
    historical_epochs = rng.random() < 0.55
    return SimulationNuisance(
        mutation_multiplier=float(np.exp(rng.uniform(np.log(0.85), np.log(1.15)))),
        recombination_multiplier=float(np.exp(rng.uniform(np.log(0.75), np.log(1.25)))),
        migration_asymmetry=float(np.exp(rng.uniform(np.log(0.25), np.log(4.0)))),
        ancestral_epoch_factors=(0.5, 1.5, 2.0)
        if historical_epochs
        else tuple(np.exp(rng.normal(0, 0.5, 3))),
        population1_epoch_factors=(0.8, 0.5)
        if historical_epochs
        else tuple(np.exp(rng.normal(0, 0.45, 2))),
        population2_epoch_factor=1.5
        if historical_epochs
        else float(np.exp(rng.normal(0, 0.4))),
        population2_contact_factor=1.0 / 3.0
        if historical_epochs
        else float(np.exp(rng.normal(0, 0.5))),
    )


def nuisance_size_changes(
    parameters: SecondaryContactParameters,
    nuisance: SimulationNuisance,
) -> list[SizeChange]:
    """Translate nuisance factors into valid hidden demographic epochs."""
    changes: list[SizeChange] = []
    for multiplier, factor in zip((1.20, 1.10, 1.05), nuisance.ancestral_epoch_factors):
        changes.append(
            SizeChange("ancestral", parameters.split_time * multiplier, parameters.ancestral_size * factor)
        )
    for fraction, factor in zip((0.80, 0.50), nuisance.population1_epoch_factors):
        time = parameters.split_time * fraction
        if time > parameters.contact_time:
            changes.append(
                SizeChange("population1", time, parameters.ancestral_size * factor)
            )
    if nuisance.population2_epoch_factor is not None:
        changes.append(
            SizeChange(
                "population2",
                parameters.contact_time * 0.45,
                parameters.population2_size * nuisance.population2_epoch_factor,
            )
        )
    if nuisance.population2_contact_factor is not None:
        changes.append(
            SizeChange(
                "population2",
                parameters.contact_time,
                parameters.population2_size * nuisance.population2_contact_factor,
            )
        )
    return changes


def simulate_block(
    parameters: SecondaryContactParameters,
    nuisance: SimulationNuisance,
    *,
    seed: int,
    sequence_length: float = 2_000_000,
    mutation_rate: float = 8.4e-8,
    recombination_rate: float = 5.39e-8,
    population1_diploids: int = 43,
    population2_diploids: int = 24,
    n_windows: int = 1,
) -> SecondaryFeatureBundle:
    """Simulate and summarize one independently seeded chromosome block."""
    demography = build_demography(
        parameters,
        size_changes=nuisance_size_changes(parameters, nuisance),
        migration_asymmetry=nuisance.migration_asymmetry,
    )
    ancestry = msprime.sim_ancestry(
        samples={"population1": population1_diploids, "population2": population2_diploids},
        ploidy=2,
        demography=demography,
        sequence_length=sequence_length,
        recombination_rate=recombination_rate * nuisance.recombination_multiplier,
        random_seed=seed,
    )
    mutated = msprime.sim_mutations(
        ancestry,
        rate=mutation_rate * nuisance.mutation_multiplier,
        model=msprime.BinaryMutationModel(),
        discrete_genome=False,
        random_seed=seed + 1,
    )
    return summarize_tree_sequence(
        mutated,
        mutated.samples(population=1),
        mutated.samples(population=2),
        n_windows=n_windows,
    )
