from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from typing import Mapping, Sequence

import numpy as np


def spectrum_from_counts(counts: Mapping[str, int], haplotypes: int) -> np.ndarray:
    """Convert a serialized unfolded SFS into a dense one-dimensional array."""
    spectrum = np.zeros(haplotypes + 1, dtype=np.float64)
    for key, value in counts.items():
        index = int(key)
        if not 0 <= index <= haplotypes:
            raise ValueError(f"SFS index {index} exceeds {haplotypes} haplotypes")
        if value < 0:
            raise ValueError("SFS counts must be nonnegative")
        spectrum[index] += value
    return spectrum


def combine_spectra(spectra: Sequence[np.ndarray]) -> np.ndarray:
    if not spectra:
        raise ValueError("at least one spectrum is required")
    arrays = [np.asarray(spectrum, dtype=np.float64) for spectrum in spectra]
    if any(array.shape != arrays[0].shape for array in arrays):
        raise ValueError("all spectra must have the same shape")
    return np.sum(arrays, axis=0)


def growth_spectrum(parameters: Sequence[float], sample_sizes: Sequence[int]):
    """Expected unfolded SFS under exponential growth.

    ``parameters`` are final/ancestral size ratio and growth duration in units
    of ``2 * N_ancestral`` generations.
    """
    import moments

    nu, duration = (float(value) for value in parameters)
    if nu <= 1 or duration <= 0:
        raise ValueError("growth ratio must exceed one and duration must be positive")
    return moments.Demographics1D.growth((nu, duration), list(sample_sizes))


def sudden_growth_spectrum(parameters: Sequence[float], sample_sizes: Sequence[int]):
    """Expected unfolded SFS after an instantaneous population-size increase."""
    import moments

    nu, duration = (float(value) for value in parameters)
    if nu <= 1 or duration <= 0:
        raise ValueError("growth ratio must exceed one and duration must be positive")
    sample_size = int(sample_sizes[0])
    spectrum = moments.Spectrum(moments.LinearSystem_1D.steady_state_1D(sample_size))
    spectrum.integrate([nu], duration)
    return spectrum


def _fit_restart(arguments: tuple) -> dict:
    spectrum, start, lower, upper, restart, maxiter, convention = arguments
    import moments

    data = moments.Spectrum(np.asarray(spectrum, dtype=np.float64))
    try:
        model = growth_spectrum if convention == "exponential" else sudden_growth_spectrum
        estimate = np.asarray(
            moments.Inference.optimize_log_lbfgsb(
                start,
                data,
                model,
                lower_bound=lower,
                upper_bound=upper,
                maxiter=maxiter,
                verbose=0,
            ),
            dtype=np.float64,
        )
        expected = model(estimate, data.sample_sizes)
        likelihood = float(moments.Inference.ll_multinom(expected, data))
        theta = float(moments.Inference.optimal_sfs_scaling(expected, data))
        if np.isfinite([*estimate, likelihood, theta]).all() and theta > 0:
            return {
                "ok": True,
                "restart": restart,
                "estimate": estimate,
                "likelihood": likelihood,
                "theta": theta,
            }
        return {"ok": False, "restart": restart, "reason": "nonfinite result"}
    except Exception as error:
        return {
            "ok": False,
            "restart": restart,
            "reason": f"{type(error).__name__}: {error}",
        }


def fit_growth_spectrum(
    spectrum: np.ndarray,
    *,
    sequence_length: float,
    mutation_rate: float = 3.36e-8,
    projection: int | None = None,
    restarts: int = 32,
    workers: int = 1,
    maxiter: int = 500,
    seed: int = 75193299,
    convention: str = "exponential",
) -> dict:
    """Fit exponential growth and convert scaled estimates to challenge units."""
    import moments

    if sequence_length <= 0 or mutation_rate <= 0:
        raise ValueError("sequence length and mutation rate must be positive")
    if restarts < 1 or workers < 1:
        raise ValueError("restarts and workers must be positive")
    if convention not in {"exponential", "sudden"}:
        raise ValueError(f"unknown growth convention: {convention}")
    data = moments.Spectrum(np.asarray(spectrum, dtype=np.float64))
    original_projection = int(data.sample_sizes[0])
    if projection is not None:
        if not 2 <= projection <= original_projection:
            raise ValueError("projection must be between 2 and the observed haplotypes")
        data = data.project([projection])
    lower = np.asarray([1.001, 1e-5], dtype=np.float64)
    upper = np.asarray([200.0, 10.0], dtype=np.float64)
    rng = np.random.default_rng(seed)
    starts = [
        np.exp(rng.uniform(np.log(lower), np.log(upper))) for _ in range(restarts)
    ]
    work = [
        (np.asarray(data), start, lower, upper, restart, maxiter, convention)
        for restart, start in enumerate(starts)
    ]
    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            results = list(executor.map(_fit_restart, work))
    else:
        results = [_fit_restart(arguments) for arguments in work]
    successful = [result for result in results if result["ok"]]
    failures = [
        {"restart": result["restart"], "reason": result["reason"]}
        for result in results
        if not result["ok"]
    ]
    if not successful:
        raise RuntimeError(f"all moments restarts failed: {failures[:3]}")
    best = max(successful, key=lambda result: result["likelihood"])
    nu, duration = best["estimate"]
    theta = best["theta"]
    ancestral_size = theta / (4 * mutation_rate * sequence_length)
    return {
        "generations": float(duration * 2 * ancestral_size),
        "post_growth_scaling": float(nu),
        "ancestral_effective_size": float(ancestral_size),
        "modern_effective_size": float(nu * ancestral_size),
        "scaled_duration_2N": float(duration),
        "theta": float(theta),
        "multinomial_log_likelihood": float(best["likelihood"]),
        "best_restart": int(best["restart"]),
        "successful_restarts": len(successful),
        "failed_restarts": failures,
        "restarts": restarts,
        "workers": workers,
        "projection": int(data.sample_sizes[0]),
        "observed_haplotypes": original_projection,
        "sequence_length": float(sequence_length),
        "mutation_rate": float(mutation_rate),
        "convention": convention,
        "bounds": {"lower": lower.tolist(), "upper": upper.tolist()},
    }
