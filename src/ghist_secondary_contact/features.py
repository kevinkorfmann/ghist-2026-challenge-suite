from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import tskit

from ghist_bottleneck.arg_features import ARG_TIME_BINS, default_time_edges
from ghist_bottleneck.features import LDAccumulator


PAIR_CLASSES = ((0, 0), (1, 1), (0, 1))
WINDOW_FEATURE_NAMES = (
    "log1p_segregating_sites_per_mb",
    "population1_pi_per_base",
    "population2_pi_per_base",
    "dxy_per_base",
    "hudson_fst",
    "private_population1_fraction",
    "private_population2_fraction",
    "shared_fraction",
    "fixed_difference_fraction",
    "mean_population1_frequency",
    "mean_population2_frequency",
    "mean_frequency_difference",
)


@dataclass(frozen=True)
class SecondaryFeatureBundle:
    windows: np.ndarray
    joint_sfs: np.ndarray
    global_summary: np.ndarray
    sequence_length: float
    n_variants: int


def population_pair_hazards(
    ts: tskit.TreeSequence,
    population1_samples: np.ndarray,
    population2_samples: np.ndarray,
    *,
    time_edges: np.ndarray | None = None,
) -> np.ndarray:
    """Return conditional coalescence hazards for three population pair classes.

    Rows represent population1-with-population1, population2-with-population2,
    and cross-population pairs. Columns are finite time intervals.
    """
    counts = population_pair_interval_counts(
        ts,
        population1_samples,
        population2_samples,
        time_edges=time_edges,
    )
    remaining = np.cumsum(counts[:, ::-1], axis=1)[:, ::-1]
    hazards = np.divide(
        counts, remaining, out=np.zeros_like(counts), where=remaining > 0
    )
    hazards = hazards[:, :-1]
    if hazards.shape != (len(PAIR_CLASSES), ARG_TIME_BINS):
        raise ValueError(f"unexpected hazard shape {hazards.shape}")
    return np.clip(hazards, 0.0, 1.0)


def population_pair_interval_counts(
    ts: tskit.TreeSequence,
    population1_samples: np.ndarray,
    population2_samples: np.ndarray,
    *,
    time_edges: np.ndarray | None = None,
) -> np.ndarray:
    """Return span-weighted coalescence counts, including the ancient tail."""
    sample_sets = [
        np.asarray(population1_samples, dtype=np.int32),
        np.asarray(population2_samples, dtype=np.int32),
    ]
    if any(len(values) == 0 for values in sample_sets):
        raise ValueError("both population sample sets must be nonempty")
    edges = (
        default_time_edges() if time_edges is None else np.asarray(time_edges, float)
    )
    if edges[0] != 0 or np.any(np.diff(edges) <= 0) or not np.all(np.isfinite(edges)):
        raise ValueError(
            "time_edges must start at zero and increase with finite values"
        )
    counts = np.asarray(
        ts.pair_coalescence_counts(
            sample_sets=sample_sets,
            indexes=list(PAIR_CLASSES),
            time_windows=np.concatenate((edges, [np.inf])),
        ),
        dtype=np.float64,
    )
    if counts.shape != (len(PAIR_CLASSES), len(edges)):
        raise ValueError(f"unexpected pair-count shape {counts.shape}")
    return counts


def transform_targets(values: np.ndarray) -> np.ndarray:
    """Map physical parameters to stable learning coordinates."""
    values = np.asarray(values, dtype=np.float64)
    if values.shape[-1] != 5:
        raise ValueError("expected five Secondary Contact parameters")
    if np.any(values <= 0):
        raise ValueError("all physical parameters must be positive")
    split = values[..., 2]
    contact = values[..., 3]
    fraction = contact / split
    if np.any((fraction <= 0) | (fraction >= 1)):
        raise ValueError("migration time must be strictly between zero and split time")
    return np.stack(
        (
            np.log(values[..., 0]),
            np.log(values[..., 1]),
            np.log(split),
            np.log(fraction) - np.log1p(-fraction),
            np.log(values[..., 4]),
        ),
        axis=-1,
    )


def inverse_transform_targets(values: np.ndarray) -> np.ndarray:
    """Map learning coordinates back to physical template parameters."""
    values = np.asarray(values, dtype=np.float64)
    if values.shape[-1] != 5:
        raise ValueError("expected five transformed parameters")
    split = np.exp(values[..., 2])
    fraction = 1.0 / (1.0 + np.exp(-values[..., 3]))
    return np.stack(
        (
            np.exp(values[..., 0]),
            np.exp(values[..., 1]),
            split,
            split * fraction,
            np.exp(values[..., 4]),
        ),
        axis=-1,
    )


def _window_row(values: np.ndarray, window_length: float) -> np.ndarray:
    sites = values[11]
    if sites <= 0:
        return np.zeros(len(WINDOW_FEATURE_NAMES), dtype=np.float64)
    fst = values[3] / max(values[2], 1e-12)
    return np.asarray(
        (
            np.log1p(sites * 1_000_000.0 / window_length),
            values[0] / window_length,
            values[1] / window_length,
            values[2] / window_length,
            fst,
            values[4] / sites,
            values[5] / sites,
            values[6] / sites,
            values[7] / sites,
            values[8] / sites,
            values[9] / sites,
            values[10] / sites,
        ),
        dtype=np.float64,
    )


def summarize_tree_sequence(
    ts: tskit.TreeSequence,
    population1_samples: np.ndarray,
    population2_samples: np.ndarray,
    *,
    n_windows: int = 50,
) -> SecondaryFeatureBundle:
    """Extract spectrum, spatial, LD, and genealogy summaries from a sequence."""
    population1_samples = np.asarray(population1_samples, dtype=np.int32)
    population2_samples = np.asarray(population2_samples, dtype=np.int32)
    if n_windows <= 0 or min(len(population1_samples), len(population2_samples)) < 2:
        raise ValueError(
            "positive window count and at least two haplotypes per population required"
        )
    sample_order = np.asarray(ts.samples(), dtype=np.int32)
    lookup = {int(node): index for index, node in enumerate(sample_order)}
    try:
        indices1 = np.asarray([lookup[int(node)] for node in population1_samples])
        indices2 = np.asarray([lookup[int(node)] for node in population2_samples])
    except KeyError as error:
        raise ValueError(
            "population set contains a node outside ts.samples()"
        ) from error
    spectrum = np.zeros((len(indices1) + 1, len(indices2) + 1), dtype=np.float64)
    # Columns store pi1, pi2, dxy, Fst numerator, four site classes,
    # frequency1, frequency2, absolute frequency difference, and site count.
    accumulators = np.zeros((n_windows, 12), dtype=np.float64)
    ld1, ld2 = LDAccumulator(seed=731), LDAccumulator(seed=947)
    window_length = ts.sequence_length / n_windows
    kept = 0
    for variant in ts.variants():
        genotype = np.asarray(variant.genotypes)
        g1, g2 = genotype[indices1], genotype[indices2]
        if np.any((g1 < 0) | (g1 > 1)) or np.any((g2 < 0) | (g2 > 1)):
            continue
        count1, count2 = int(g1.sum()), int(g2.sum())
        if count1 + count2 in (0, len(g1) + len(g2)):
            continue
        spectrum[count1, count2] += 1
        p1, p2 = count1 / len(g1), count2 / len(g2)
        pi1 = 2 * p1 * (1 - p1) * len(g1) / (len(g1) - 1)
        pi2 = 2 * p2 * (1 - p2) * len(g2) / (len(g2) - 1)
        dxy = p1 * (1 - p2) + p2 * (1 - p1)
        row = min(int(variant.site.position / window_length), n_windows - 1)
        values = accumulators[row]
        values[0] += pi1
        values[1] += pi2
        values[2] += dxy
        values[3] += dxy - 0.5 * (pi1 + pi2)
        values[4] += int(0 < count1 < len(g1) and count2 in (0, len(g2)))
        values[5] += int(0 < count2 < len(g2) and count1 in (0, len(g1)))
        values[6] += int(0 < count1 < len(g1) and 0 < count2 < len(g2))
        values[7] += int(
            (count1 == 0 and count2 == len(g2)) or (count2 == 0 and count1 == len(g1))
        )
        values[8] += p1
        values[9] += p2
        values[10] += abs(p1 - p2)
        values[11] += 1
        position = float(variant.site.position)
        paired1 = g1[: (len(g1) // 2) * 2].reshape(-1, 2).sum(axis=1)
        paired2 = g2[: (len(g2) // 2) * 2].reshape(-1, 2).sum(axis=1)
        ld1.add(position, paired1)
        ld2.add(position, paired2)
        kept += 1
    windows = np.vstack([_window_row(values, window_length) for values in accumulators])
    total = accumulators.sum(axis=0)
    global_spatial = _window_row(total, ts.sequence_length)
    hazards = population_pair_hazards(ts, population1_samples, population2_samples)
    global_summary = np.concatenate(
        (global_spatial, ld1.result(), ld2.result(), hazards.reshape(-1))
    )
    spectrum = np.log1p(spectrum)
    spectrum /= np.linalg.norm(spectrum) + 1e-12
    return SecondaryFeatureBundle(
        windows=windows,
        joint_sfs=spectrum,
        global_summary=global_summary,
        sequence_length=float(ts.sequence_length),
        n_variants=kept,
    )
