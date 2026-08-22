from __future__ import annotations

from dataclasses import asdict

import numpy as np
from scipy.linalg import expm
from scipy.optimize import differential_evolution, minimize

from ghist_bottleneck.arg_features import default_time_edges

from .demography import SecondaryContactParameters


PAIR_ORDER = (0, 2, 1)  # mainland/mainland, island/island, mainland/island


def _validate_edges(time_edges: np.ndarray) -> np.ndarray:
    edges = np.asarray(time_edges, dtype=np.float64)
    if edges.ndim != 1 or edges[0] != 0 or np.any(np.diff(edges) <= 0):
        raise ValueError("time edges must be one-dimensional, start at zero, and increase")
    if not np.all(np.isfinite(edges)):
        raise ValueError("finite time edges are required")
    return edges


def secondary_contact_interval_probabilities(
    parameters: SecondaryContactParameters,
    *,
    time_edges: np.ndarray | None = None,
) -> np.ndarray:
    """Exact pairwise coalescence masses for a symmetric secondary-contact model.

    The three transient states are mainland/mainland, cross-population, and
    island/island. The returned final column is the tail beyond the last finite
    time edge, matching ``tskit.pair_coalescence_counts``.
    """
    parameters.validate()
    edges = _validate_edges(default_time_edges() if time_edges is None else time_edges)
    n1 = parameters.population1_size
    n2 = parameters.population2_size
    ancestral = parameters.ancestral_size
    migration = parameters.migration_rate
    contact = parameters.contact_time
    split = parameters.split_time
    c1, c2, ca = 1.0 / (2.0 * n1), 1.0 / (2.0 * n2), 1.0 / (2.0 * ancestral)
    contact_generator = np.asarray(
        (
            (-(2.0 * migration + c1), 2.0 * migration, 0.0),
            (migration, -2.0 * migration, migration),
            (0.0, 2.0 * migration, -(2.0 * migration + c2)),
        ),
        dtype=np.float64,
    )
    isolation_generator = np.diag((-c1, 0.0, -c2))
    contact_transition = expm(contact_generator * contact)
    isolation_transition = expm(isolation_generator * (split - contact))
    survival = np.empty((3, len(edges)), dtype=np.float64)
    identity = np.eye(3, dtype=np.float64)
    for column, time in enumerate(edges):
        if time <= contact:
            transition = expm(contact_generator * time)
            survival[:, column] = transition.sum(axis=1)
        elif time <= split:
            transition = contact_transition @ expm(
                isolation_generator * (time - contact)
            )
            survival[:, column] = transition.sum(axis=1)
        else:
            at_split = (contact_transition @ isolation_transition).sum(axis=1)
            survival[:, column] = at_split * np.exp(-ca * (time - split))
    survival[:, 0] = identity.sum(axis=1)
    finite_mass = np.clip(survival[:, :-1] - survival[:, 1:], 0.0, 1.0)
    tail = np.clip(survival[:, -1:], 0.0, 1.0)
    probabilities = np.concatenate((finite_mass, tail), axis=1)[list(PAIR_ORDER)]
    probabilities /= probabilities.sum(axis=1, keepdims=True)
    return probabilities


def secondary_contact_pair_hazards(
    parameters: SecondaryContactParameters,
    *,
    time_edges: np.ndarray | None = None,
) -> np.ndarray:
    probabilities = secondary_contact_interval_probabilities(
        parameters,
        time_edges=time_edges,
    )
    remaining = np.cumsum(probabilities[:, ::-1], axis=1)[:, ::-1]
    return np.divide(
        probabilities[:, :-1],
        remaining[:, :-1],
        out=np.zeros_like(probabilities[:, :-1]),
        where=remaining[:, :-1] > 0,
    )


def _from_coordinates(coordinates: np.ndarray) -> SecondaryContactParameters:
    split = float(np.exp(coordinates[3]))
    fraction = float(1.0 / (1.0 + np.exp(-coordinates[4])))
    return SecondaryContactParameters(
        ancestral_size=float(np.exp(coordinates[2])),
        population1_size=float(np.exp(coordinates[0])),
        population2_size=float(np.exp(coordinates[1])),
        split_time=split,
        contact_time=split * fraction,
        migration_rate=float(np.exp(coordinates[5])),
    )


def _to_coordinates(parameters: SecondaryContactParameters) -> np.ndarray:
    fraction = parameters.contact_time / parameters.split_time
    return np.asarray(
        (
            np.log(parameters.population1_size),
            np.log(parameters.population2_size),
            np.log(parameters.ancestral_size),
            np.log(parameters.split_time),
            np.log(fraction) - np.log1p(-fraction),
            np.log(parameters.migration_rate),
        ),
        dtype=np.float64,
    )


def fit_secondary_contact_interval_counts(
    interval_counts: np.ndarray,
    *,
    initial: SecondaryContactParameters,
    global_search: bool = True,
    seed: int = 20260822,
) -> dict:
    """Fit physical parameters to span-weighted pairwise coalescence counts."""
    counts = np.asarray(interval_counts, dtype=np.float64)
    expected_shape = (3, len(default_time_edges()))
    if counts.shape != expected_shape or np.any(counts < 0):
        raise ValueError(f"interval counts must have shape {expected_shape} and be nonnegative")
    totals = counts.sum(axis=1, keepdims=True)
    if np.any(totals <= 0):
        raise ValueError("every pair class must have positive total weight")
    empirical = counts / totals
    bounds = (
        (np.log(1_000.0), np.log(300_000.0)),
        (np.log(1_000.0), np.log(300_000.0)),
        (np.log(1_000.0), np.log(300_000.0)),
        (np.log(100.0), np.log(200_000.0)),
        (-8.0, 8.0),
        (np.log(1e-8), np.log(0.05)),
    )

    def objective(coordinates: np.ndarray) -> float:
        probabilities = secondary_contact_interval_probabilities(
            _from_coordinates(coordinates)
        )
        return float(-np.sum(empirical * np.log(np.clip(probabilities, 1e-15, 1.0))))

    candidates = [(_to_coordinates(initial), "mechanistic_initial")]
    if global_search:
        broad = differential_evolution(
            objective,
            bounds=bounds,
            seed=seed,
            popsize=10,
            maxiter=100,
            polish=False,
            workers=1,
            updating="immediate",
        )
        candidates.append((broad.x, "differential_evolution"))
    fits = []
    for start, label in candidates:
        fit = minimize(objective, start, method="L-BFGS-B", bounds=bounds)
        fits.append((float(fit.fun), fit, label))
    loss, fit, source = min(fits, key=lambda item: item[0])
    parameters = _from_coordinates(fit.x)
    predicted = secondary_contact_interval_probabilities(parameters)
    hazards = secondary_contact_pair_hazards(parameters)
    return {
        **asdict(parameters),
        "negative_equal_pair_class_log_likelihood": loss,
        "converged": bool(fit.success),
        "optimizer_message": str(fit.message),
        "initialization": source,
        "empirical_interval_probabilities": empirical.tolist(),
        "predicted_interval_probabilities": predicted.tolist(),
        "predicted_hazards": hazards.tolist(),
    }
