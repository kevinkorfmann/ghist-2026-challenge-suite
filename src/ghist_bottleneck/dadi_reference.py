from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from .sfs_baseline import projected_sfs_from_vcf, spectrum_exposure
from .optimization import DEFAULT_STRATEGY, STRATEGIES, LOWER, UPPER, log_starts, run_optimizer
from ghist_inference.validation import (
    validate_spectrum, validate_fit_settings, validate_expected_spectrum, optimizer_diagnostics,
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
    initial: tuple[float, float] | None = None,
    optimizer_strategy: str = DEFAULT_STRATEGY,
    initial_log_starts: np.ndarray | None = None,
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
    def evaluate(parameters):
        theta = float(dadi.Inference.optimal_sfs_scaling(model_function(parameters, data.sample_sizes, grids), data))
        if not np.isfinite(theta) or theta <= 0:
            raise ValueError("Nonfinite or nonpositive theta")
        return {"theta": theta}

    starts = log_starts(restarts, seed, initial=initial, starts=initial_log_starts)
    optimization = run_optimizer(objective, engine="dadi", starts=starts, strategy=optimizer_strategy, evaluate=evaluate)
    best = optimization["best"]
    estimate = np.asarray(best["parameters"])
    likelihood, theta, best_restart = -best["objective"], best["theta"], best["restart"]
    optimizer = best["optimizer"]
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
        "restarts": restarts,
        "engine": "dadi",
        "optimizer_contract": optimization["contract"],
        "sample_haplotypes": int(data.sample_sizes[0]),
        "optimizer_strategy": optimizer_strategy,
        "best_optimizer_arm": best["arm"],
        "best_restart_converged": optimizer["converged"],
        "successful_restarts": len(optimization["origin_best"]),
        "converged_restarts": sum(row["optimizer"]["converged"] for row in optimization["origin_best"]),
        "failed_restarts": optimization["failed_origins"],
        "grid_points": grids,
        "optimizer": optimizer,
        "optimizer_neighborhood": best["neighborhood"],
        "optimizer_attempts": optimization["attempts"],
        "optimizer_arm_summary": optimization["arm_summary"],
        "optimizer_seed": seed,
        "optimizer_initialization": "explicit_log_starts" if initial_log_starts is not None else "seeded_random_with_optional_first_initial",
        "optimizer_origin_best": [{"restart": row["restart"], "arm": row["arm"], "objective": row["objective"],
                                   "converged": row["optimizer"]["converged"]} for row in optimization["origin_best"]],
        "optimizer_diagnostics": optimizer_diagnostics(
            estimate, LOWER, UPPER, [-row["objective"] for row in optimization["origin_best"]], likelihood
        ),
        "folded": folded,
        "mutation_rate": mutation_rate,
        "sequence_length": sequence_length,
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
    optimizer_strategy: str = DEFAULT_STRATEGY,
) -> dict:
    with Path(vcf).open("rb") as handle:
        input_sha256 = hashlib.file_digest(handle, "sha256").hexdigest()
    spectrum, metadata = projected_sfs_from_vcf(vcf, polarized=not folded)
    exposure = spectrum_exposure(metadata, sequence_length, callable_sequence_length)
    result = fit_two_epoch_dadi_spectrum(
        spectrum,
        restarts=restarts,
        folded=folded,
        seed=seed,
        mutation_rate=mutation_rate,
        sequence_length=exposure,
        optimizer_strategy=optimizer_strategy,
    )
    with Path(vcf).open("rb") as handle:
        if hashlib.file_digest(handle, "sha256").hexdigest() != input_sha256:
            raise ValueError("VCF changed during spectrum extraction or fitting; regenerate the fit")
    return {**result, "polarized": not folded, "input_vcf_sha256": input_sha256,
            "sfs": metadata, "sequence_length": sequence_length,
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
    parser.add_argument("--optimizer-strategy", choices=STRATEGIES, default=DEFAULT_STRATEGY)
    args = parser.parse_args()
    result = fit_two_epoch_dadi(
        args.vcf,
        restarts=args.restarts,
        folded=args.folded,
        seed=args.seed,
        sequence_length=args.sequence_length,
        mutation_rate=args.mutation_rate,
        callable_sequence_length=args.callable_sequence_length,
        optimizer_strategy=args.optimizer_strategy,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
