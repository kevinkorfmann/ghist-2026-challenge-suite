from __future__ import annotations

import numpy as np
import tskit


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
