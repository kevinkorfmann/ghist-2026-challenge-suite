from __future__ import annotations

import numpy as np
import tskit
from scipy.special import ndtr


ARG_TIME_BINS = 32


def default_time_edges(n_bins: int = ARG_TIME_BINS) -> np.ndarray:
    """Edges covering recent to moderately ancient demographic events."""
    return np.concatenate(([0.0], np.geomspace(25.0, 200_000.0, n_bins)))


def coalescence_hazard_features(
    ts: tskit.TreeSequence,
    *,
    samples: np.ndarray | None = None,
    time_edges: np.ndarray | None = None,
) -> np.ndarray:
    """Conditional interval coalescence probabilities from a tree sequence.

    This follows the coalescent-based time-stratification principle of Deng,
    Pritchard, and Spence (2026): counts in each interval are divided by the
    number of lineages/pairs that have survived without coalescing to its start.
    """
    if samples is None:
        samples = ts.samples()
    samples = np.asarray(samples, dtype=np.int32)
    if len(samples) < 2:
        return np.zeros(ARG_TIME_BINS, dtype=np.float64)
    edges = default_time_edges() if time_edges is None else np.asarray(time_edges, dtype=float)
    if edges[0] != 0 or np.any(np.diff(edges) <= 0) or not np.all(np.isfinite(edges)):
        raise ValueError("time_edges must start at zero and increase with finite values")
    tskit_edges = np.concatenate((edges, [np.inf]))
    counts = np.asarray(
        ts.pair_coalescence_counts(sample_sets=[samples], time_windows=tskit_edges),
        dtype=np.float64,
    ).reshape(-1)
    remaining = np.cumsum(counts[::-1])[::-1]
    hazard = np.divide(counts, remaining, out=np.zeros_like(counts), where=remaining > 0)
    # The infinity tail is required by tskit but is not a finite-time feature.
    hazard = hazard[:-1]
    # Fixed length protects serialized datasets if custom edges are accidentally used.
    if len(hazard) != ARG_TIME_BINS:
        raise ValueError(f"expected {ARG_TIME_BINS} hazard bins, got {len(hazard)}")
    return hazard


def tmrca_hazard_features(
    log_tmrca: np.ndarray,
    *,
    time_edges: np.ndarray | None = None,
) -> np.ndarray:
    """Convert posterior log-TMRCA samples into conditional time hazards."""
    values = np.exp(np.asarray(log_tmrca, dtype=np.float64).reshape(-1))
    values = values[np.isfinite(values) & (values >= 0)]
    if not len(values):
        return np.zeros(ARG_TIME_BINS, dtype=np.float64)
    edges = default_time_edges() if time_edges is None else np.asarray(time_edges, dtype=float)
    if edges[0] != 0 or np.any(np.diff(edges) <= 0) or not np.all(np.isfinite(edges)):
        raise ValueError("time_edges must start at zero and increase with finite values")
    counts, _ = np.histogram(values, bins=np.concatenate((edges, [np.inf])))
    counts = counts.astype(np.float64)
    remaining = np.cumsum(counts[::-1])[::-1]
    hazard = np.divide(counts, remaining, out=np.zeros_like(counts), where=remaining > 0)
    hazard = hazard[:-1]
    if len(hazard) != ARG_TIME_BINS:
        raise ValueError(f"expected {ARG_TIME_BINS} hazard bins, got {len(hazard)}")
    return hazard


def gaussian_tmrca_hazard_features(
    log_mean: np.ndarray,
    log_variance: np.ndarray,
    *,
    time_edges: np.ndarray | None = None,
) -> np.ndarray:
    """Integrate Gaussian log-TMRCA posteriors into conditional hazards.

    Each pair/window prediction defines a Gaussian distribution in log time.
    Interval probabilities are evaluated analytically and averaged before the
    conditional-hazard transform, avoiding posterior Monte Carlo noise.
    """
    means = np.asarray(log_mean, dtype=np.float64)
    variances = np.asarray(log_variance, dtype=np.float64)
    if means.shape != variances.shape:
        raise ValueError(
            f"log_mean and log_variance must have the same shape, got "
            f"{means.shape} and {variances.shape}"
        )
    means = means.reshape(-1)
    variances = variances.reshape(-1)
    keep = np.isfinite(means) & np.isfinite(variances) & (variances >= 0)
    means, variances = means[keep], variances[keep]
    if not len(means):
        return np.zeros(ARG_TIME_BINS, dtype=np.float64)

    edges = default_time_edges() if time_edges is None else np.asarray(time_edges, dtype=float)
    if edges[0] != 0 or np.any(np.diff(edges) <= 0) or not np.all(np.isfinite(edges)):
        raise ValueError("time_edges must start at zero and increase with finite values")

    # CDF values at the finite positive edges. The lower support boundary is
    # exactly zero (log boundary -infinity), and the final tail ends at infinity.
    log_bounds = np.log(edges[1:])
    sigma = np.sqrt(variances)
    cdf = np.empty((len(means), len(log_bounds)), dtype=np.float64)
    uncertain = sigma > 1e-8
    if np.any(uncertain):
        z = (log_bounds[None, :] - means[uncertain, None]) / sigma[uncertain, None]
        cdf[uncertain] = ndtr(z)
    if np.any(~uncertain):
        cdf[~uncertain] = means[~uncertain, None] <= log_bounds[None, :]

    interval_mass = np.diff(
        np.concatenate(
            (
                np.zeros((len(means), 1), dtype=np.float64),
                cdf,
                np.ones((len(means), 1), dtype=np.float64),
            ),
            axis=1,
        ),
        axis=1,
    ).mean(axis=0)
    interval_mass = np.clip(interval_mass, 0.0, 1.0)
    remaining = np.cumsum(interval_mass[::-1])[::-1]
    hazard = np.divide(
        interval_mass,
        remaining,
        out=np.zeros_like(interval_mass),
        where=remaining > 0,
    )[:-1]
    if len(hazard) != ARG_TIME_BINS:
        raise ValueError(f"expected {ARG_TIME_BINS} hazard bins, got {len(hazard)}")
    return np.clip(hazard, 0.0, 1.0)
