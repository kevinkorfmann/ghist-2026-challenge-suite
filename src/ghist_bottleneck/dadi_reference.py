from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.optimize import fmin_l_bfgs_b

from .sfs_baseline import projected_sfs_from_vcf, spectrum_exposure
from ghist_inference.validation import (
    validate_spectrum, validate_fit_settings, validate_expected_spectrum,
)


def physical_parameters(
    fraction: float,
    scaled_time: float,
    theta: float,
    *,
    mutation_rate: float,
    sequence_length: float,
) -> tuple[float, float, float]:
    ancestral_ne = theta / (4 * mutation_rate * sequence_length)
    generations = scaled_time * 2 * ancestral_ne
    return float(generations), float(fraction), float(ancestral_ne)


def fit_two_epoch_dadi_spectrum(
    spectrum: np.ndarray,
    *,
    mutation_rate: float = 3.42e-8,
    sequence_length: float = 100_000_000,
    restarts: int = 64,
    seed: int = 8241,
    folded: bool = False,
) -> dict:
    """Fit dadi's diffusion solution with independent log-space restarts."""
    import dadi

    spectrum = validate_spectrum(spectrum, 1)
    validate_fit_settings(
        sequence_length=sequence_length, mutation_rate=mutation_rate,
        restarts=restarts, workers=1, maxiter=500, shape=spectrum.shape,
    )
    data = dadi.Spectrum(spectrum)
    if folded:
        data = data.fold()
    sample_size = int(data.sample_sizes[0])
    grids = [sample_size + offset for offset in (20, 30, 40)]
    extrapolated = dadi.Numerics.make_extrap_log_func(dadi.Demographics1D.two_epoch)

    def model_function(params, ns, pts):
        expected = extrapolated(params, ns, pts)
        expected = expected.fold() if folded else expected
        return validate_expected_spectrum(expected, data)

    def objective(log_parameters):
        expected = model_function(np.exp(log_parameters), data.sample_sizes, grids)
        likelihood = float(dadi.Inference.ll_multinom(expected, data))
        if not np.isfinite(likelihood):
            raise ValueError("Nonfinite diffusion likelihood")
        return -likelihood
    lower = np.array([0.02, 1e-4])
    upper = np.array([0.999, 5.0])
    rng = np.random.default_rng(seed)
    fits = []
    failures = []
    for restart in range(restarts):
        start = np.exp(rng.uniform(np.log(lower), np.log(upper)))
        try:
            # dadi 2.4's wrapper passes removed SciPy `iprint`. Use the same
            # public optimizer and log objective without that obsolete option.
            optimum, _, optimizer = fmin_l_bfgs_b(
                objective, np.log(start),
                bounds=list(zip(np.log(lower), np.log(upper))),
                epsilon=1e-3, pgtol=1e-5, maxfun=500, approx_grad=True,
            )
            estimate = np.exp(optimum)
            expected = model_function(estimate, data.sample_sizes, grids)
            likelihood = float(dadi.Inference.ll_multinom(expected, data))
            theta = float(dadi.Inference.optimal_sfs_scaling(expected, data))
            if np.isfinite([*estimate, likelihood, theta]).all() and theta > 0:
                fits.append((likelihood, estimate, theta, restart, optimizer))
            else:
                failures.append({"restart": restart, "reason": "nonfinite result"})
        except Exception as error:
            failures.append(
                {"restart": restart, "reason": f"{type(error).__name__}: {error}"}
            )
    if not fits:
        raise RuntimeError(f"all dadi restarts failed: {failures[:3]}")
    likelihood, estimate, theta, best_restart, optimizer = max(fits, key=lambda item: item[0])
    fraction, scaled_time = estimate
    generations, fraction, ancestral_ne = physical_parameters(
        fraction,
        scaled_time,
        theta,
        mutation_rate=mutation_rate,
        sequence_length=sequence_length,
    )
    return {
        "generations": generations,
        "post_decline_fraction": fraction,
        "ancestral_effective_size": ancestral_ne,
        "scaled_time_2N": float(scaled_time),
        "theta": theta,
        "multinomial_log_likelihood": likelihood,
        "best_restart": int(best_restart),
        "successful_restarts": len(fits),
        "failed_restarts": failures,
        "grid_points": grids,
        "optimizer": {"warnflag": int(optimizer["warnflag"]),
                      "message": str(optimizer["task"]),
                      "evaluations": int(optimizer["funcalls"])},
        "folded": folded,
    }


def fit_two_epoch_dadi(
    vcf: str | Path,
    *,
    restarts: int = 64,
    folded: bool = False,
    seed: int = 8241,
    mutation_rate: float = 3.42e-8,
    sequence_length: float = 100_000_000,
    callable_sequence_length: float | None = None,
) -> dict:
    spectrum, metadata = projected_sfs_from_vcf(vcf, polarized=not folded)
    exposure = spectrum_exposure(metadata, sequence_length, callable_sequence_length)
    result = fit_two_epoch_dadi_spectrum(
        spectrum,
        restarts=restarts,
        folded=folded,
        seed=seed,
        mutation_rate=mutation_rate,
        sequence_length=exposure,
    )
    return {**result, "sfs": metadata, "sequence_length": sequence_length,
            "callable_sequence_length": exposure}


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--restarts", type=int, default=64)
    parser.add_argument("--seed", type=int, default=8241)
    parser.add_argument("--folded", action="store_true")
    parser.add_argument("--sequence-length", type=float, default=100_000_000)
    parser.add_argument("--mutation-rate", type=float, default=3.42e-8)
    parser.add_argument("--callable-sequence-length", type=float)
    args = parser.parse_args()
    result = fit_two_epoch_dadi(
        args.vcf,
        restarts=args.restarts,
        folded=args.folded,
        seed=args.seed,
        sequence_length=args.sequence_length,
        mutation_rate=args.mutation_rate,
        callable_sequence_length=args.callable_sequence_length,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
