from __future__ import annotations

import numpy as np


def validate_expected_spectrum(expected, data):
    """Keep the likelihood domain fixed when a numerical integration fails.

    moments/dadi masked-array likelihoods can silently omit negative, nonfinite,
    zero-probability, or additionally masked model bins. A finite score is not
    sufficient evidence that every observed bin contributed to the likelihood.
    """
    observed = np.asarray(data, dtype=float)
    predicted = np.asarray(expected, dtype=float)
    if predicted.shape != observed.shape:
        raise ValueError("Expected and observed spectra have different shapes")
    active = ~np.ma.getmaskarray(data)
    invalid = active & (
        np.ma.getmaskarray(expected)
        | ~np.isfinite(predicted)
        | (predicted < 0)
        | ((predicted == 0) & (observed > 0))
    )
    if np.any(invalid):
        raise ValueError("Expected spectrum would exclude observed likelihood bins")
    return expected


def validate_spectrum(spectrum, dimensions: int) -> np.ndarray:
    values = np.asarray(spectrum, dtype=float)
    if (
        values.ndim != dimensions
        or any(size < 3 for size in values.shape)
        or not np.isfinite(values).all()
        or np.any(values < 0)
    ):
        raise ValueError(
            f"Expected a finite nonnegative {dimensions}D spectrum with at least two haplotypes per population"
        )
    if values.ravel()[1:-1].sum() <= 0:
        raise ValueError("Spectrum has no polymorphic sites")
    return values


def validate_fit_settings(
    *,
    sequence_length,
    mutation_rate,
    restarts,
    workers,
    maxiter,
    shape,
    projection=None,
    integration_dt_fac=0.02,
):
    if (
        not np.isfinite([sequence_length, mutation_rate, integration_dt_fac]).all()
        or min(sequence_length, mutation_rate, integration_dt_fac) <= 0
    ):
        raise ValueError(
            "Sequence length, mutation rate, and integration step must be finite and positive"
        )
    for name, value in (
        ("restarts", restarts),
        ("workers", workers),
        ("maxiter", maxiter),
    ):
        if not isinstance(value, (int, np.integer)) or value < 1:
            raise ValueError(f"{name} must be a positive integer")
    if projection is not None:
        if len(projection) != len(shape) or any(
            not isinstance(n, (int, np.integer)) or not 2 <= n < size
            for n, size in zip(projection, shape)
        ):
            raise ValueError(
                "Projection must contain one integer per population, between 2 and its observed haplotypes"
            )


def optimizer_diagnostics(estimate, lower, upper, likelihoods, best_likelihood):
    estimate, lower, upper = (
        np.asarray(x, dtype=float) for x in (estimate, lower, upper)
    )
    near = (estimate <= lower * 1.01) | (estimate >= upper / 1.01)
    support = int(
        np.count_nonzero(np.abs(np.asarray(likelihoods) - best_likelihood) <= 0.1)
    )
    return {
        "near_boundary_parameter_indices": np.flatnonzero(near).tolist(),
        "restarts_within_0_1_log_likelihood": support,
        "numerically_supported": bool(not near.any() and support >= 2),
        "scope": "Restart agreement and boundary diagnostic; not proof of model adequacy or optimizer convergence",
    }
