from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from .sfs_baseline import projected_sfs_from_vcf


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

    data = dadi.Spectrum(np.asarray(spectrum, dtype=float))
    if folded:
        data = data.fold()
    sample_size = int(data.sample_sizes[0])
    grids = [sample_size + offset for offset in (20, 30, 40)]
    model_function = dadi.Numerics.make_extrap_log_func(dadi.Demographics1D.two_epoch)
    lower = np.array([0.02, 1e-4])
    upper = np.array([0.999, 5.0])
    rng = np.random.default_rng(seed)
    fits = []
    failures = []
    for restart in range(restarts):
        start = np.exp(rng.uniform(np.log(lower), np.log(upper)))
        try:
            estimate = np.asarray(
                dadi.Inference.optimize_log_lbfgsb(
                    start,
                    data,
                    model_function,
                    grids,
                    lower_bound=lower,
                    upper_bound=upper,
                    maxiter=500,
                    verbose=0,
                ),
                dtype=float,
            )
            expected = model_function(estimate, data.sample_sizes, grids)
            likelihood = float(dadi.Inference.ll_multinom(expected, data))
            theta = float(dadi.Inference.optimal_sfs_scaling(expected, data))
            if np.isfinite([*estimate, likelihood, theta]).all():
                fits.append((likelihood, estimate, theta, restart))
            else:
                failures.append({"restart": restart, "reason": "nonfinite result"})
        except Exception as error:
            failures.append(
                {"restart": restart, "reason": f"{type(error).__name__}: {error}"}
            )
    if not fits:
        raise RuntimeError(f"all dadi restarts failed: {failures[:3]}")
    likelihood, estimate, theta, best_restart = max(fits, key=lambda item: item[0])
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
        "folded": folded,
    }


def fit_two_epoch_dadi(
    vcf: str | Path,
    *,
    restarts: int = 64,
    folded: bool = False,
    seed: int = 8241,
) -> dict:
    spectrum, metadata = projected_sfs_from_vcf(vcf, polarized=not folded)
    result = fit_two_epoch_dadi_spectrum(
        spectrum,
        restarts=restarts,
        folded=folded,
        seed=seed,
    )
    return {**result, "sfs": metadata}


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--restarts", type=int, default=64)
    parser.add_argument("--seed", type=int, default=8241)
    parser.add_argument("--folded", action="store_true")
    args = parser.parse_args()
    result = fit_two_epoch_dadi(
        args.vcf,
        restarts=args.restarts,
        folded=args.folded,
        seed=args.seed,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
