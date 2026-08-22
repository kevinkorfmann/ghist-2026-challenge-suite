from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import pysam
from scipy.stats import hypergeom


def secondary_contact_spectrum(
    parameters: Sequence[float],
    sample_sizes: Sequence[int],
    *,
    folded: bool = False,
):
    """Expected joint SFS for isolation followed by symmetric migration.

    Parameters are population-size ratios, isolation duration, contact duration,
    and scaled migration: ``nu1, nu2, Tiso, Tcontact, M``. Times use units of
    ``2 * Nref`` and ``M = 2 * Nref * m``.
    """
    import moments

    nu1, nu2, isolation_duration, contact_duration, migration = parameters
    if isolation_duration <= 0 or contact_duration <= 0:
        raise ValueError("epoch durations must be positive")
    n1, n2 = (int(value) for value in sample_sizes)
    spectrum = moments.Spectrum(
        moments.LinearSystem_1D.steady_state_1D(n1 + n2)
    )
    spectrum = moments.Manips.split_1D_to_2D(spectrum, n1, n2)
    spectrum.integrate([nu1, nu2], isolation_duration)
    migration_matrix = np.asarray([[0.0, migration], [migration, 0.0]])
    spectrum.integrate(
        [nu1, nu2],
        contact_duration,
        m=migration_matrix,
    )
    return spectrum.fold() if folded else spectrum


def _fit_restart(arguments: tuple) -> dict:
    spectrum, start, lower, upper, folded, restart, maxiter = arguments
    import moments

    data = moments.Spectrum(np.asarray(spectrum, dtype=float))
    if folded:
        data = data.fold()

    def model(parameters, sample_sizes):
        return secondary_contact_spectrum(parameters, sample_sizes, folded=folded)

    try:
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
            dtype=float,
        )
        expected = model(estimate, data.sample_sizes)
        likelihood = float(moments.Inference.ll_multinom(expected, data))
        theta = float(moments.Inference.optimal_sfs_scaling(expected, data))
        if np.isfinite([*estimate, likelihood, theta]).all():
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


def fit_secondary_contact_spectrum(
    spectrum: np.ndarray,
    *,
    mutation_rate: float = 8.4e-8,
    sequence_length: float = 100_000_000,
    restarts: int = 32,
    workers: int = 1,
    projection: tuple[int, int] | None = None,
    maxiter: int = 750,
    seed: int = 74031,
    folded: bool = False,
) -> dict:
    """Fit the five-parameter joint-SFS reference with independent restarts."""
    fitting_spectrum = np.asarray(spectrum, dtype=float)
    if projection is not None:
        import moments

        fitting_spectrum = np.asarray(
            moments.Spectrum(fitting_spectrum).project(list(projection)), dtype=float
        )
    lower = np.asarray([0.01, 0.01, 1e-5, 1e-5, 1e-6])
    upper = np.asarray([100.0, 100.0, 10.0, 10.0, 100.0])
    rng = np.random.default_rng(seed)
    starts = [np.exp(rng.uniform(np.log(lower), np.log(upper))) for _ in range(restarts)]
    work = [
        (fitting_spectrum, start, lower, upper, folded, restart, maxiter)
        for restart, start in enumerate(starts)
    ]
    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            results = list(executor.map(_fit_restart, work))
    else:
        results = [_fit_restart(arguments) for arguments in work]
    fits = [
        (result["likelihood"], result["estimate"], result["theta"], result["restart"])
        for result in results
        if result["ok"]
    ]
    failures = [
        {"restart": result["restart"], "reason": result["reason"]}
        for result in results
        if not result["ok"]
    ]
    if not fits:
        raise RuntimeError(f"all moments restarts failed: {failures[:3]}")
    likelihood, estimate, theta, best_restart = max(fits, key=lambda item: item[0])
    nu1, nu2, isolation_duration, contact_duration, scaled_migration = estimate
    ancestral_size = theta / (4 * mutation_rate * sequence_length)
    contact_time = contact_duration * 2 * ancestral_size
    split_time = (isolation_duration + contact_duration) * 2 * ancestral_size
    migration_rate = scaled_migration / (2 * ancestral_size)
    return {
        "generations_since_split": float(split_time),
        "generations_since_contact": float(contact_time),
        "population1_effective_size": float(nu1 * ancestral_size),
        "population2_effective_size": float(nu2 * ancestral_size),
        "ancestral_effective_size": float(ancestral_size),
        "migration_rate": float(migration_rate),
        "scaled_parameters": {
            "population1_size_ratio": float(nu1),
            "population2_size_ratio": float(nu2),
            "isolation_duration_2N": float(isolation_duration),
            "contact_duration_2N": float(contact_duration),
            "migration_2N": float(scaled_migration),
        },
        "theta": theta,
        "multinomial_log_likelihood": likelihood,
        "best_restart": int(best_restart),
        "successful_restarts": len(fits),
        "failed_restarts": failures,
        "restarts": restarts,
        "workers": workers,
        "projection": None if projection is None else list(projection),
        "maxiter": maxiter,
        "folded": folded,
    }


def joint_spectrum_from_vcf(
    path: str | Path,
    populations: Mapping[str, Sequence[str]],
    *,
    projection: tuple[int, int] | None = None,
    polarized: bool = True,
) -> tuple[np.ndarray, dict]:
    """Extract an exact or hypergeometrically projected two-population SFS."""
    if len(populations) != 2:
        raise ValueError("exactly two populations are required")
    names = list(populations)
    sample_groups = [list(populations[name]) for name in names]
    vcf = pysam.VariantFile(str(path))
    header_samples = set(vcf.header.samples)
    requested = [sample for group in sample_groups for sample in group]
    if len(requested) != len(set(requested)):
        raise ValueError("a sample cannot belong to both populations")
    missing = sorted(set(requested) - header_samples)
    if missing:
        raise ValueError(f"samples absent from VCF: {missing[:5]}")
    maximum = tuple(2 * len(group) for group in sample_groups)
    if projection is None:
        projection = maximum
    if any(value < 2 or value > limit for value, limit in zip(projection, maximum)):
        raise ValueError(f"invalid projection {projection}; maxima are {maximum}")
    result = np.zeros((projection[0] + 1, projection[1] + 1), dtype=float)
    usable = projected = ancestral_reference = ancestral_alternate = 0
    for record in vcf:
        if len(record.alleles or ()) != 2 or any(len(a) != 1 for a in record.alleles):
            continue
        derived = []
        called = []
        for group in sample_groups:
            alleles = []
            for sample in group:
                genotype = record.samples[sample].get("GT")
                if genotype is not None:
                    alleles.extend(a for a in genotype if a is not None and a >= 0)
            called.append(len(alleles))
            derived.append(sum(a == 1 for a in alleles))
        if polarized:
            ancestral = record.info.get("AA")
            if isinstance(ancestral, tuple):
                ancestral = ancestral[0]
            if ancestral == record.ref:
                ancestral_reference += 1
            elif ancestral == record.alts[0]:
                ancestral_alternate += 1
                derived = [n - count for n, count in zip(called, derived)]
            else:
                raise ValueError(
                    f"record {record.contig}:{record.pos} has invalid ancestral allele"
                )
        if all(n >= target for n, target in zip(called, projection)):
            weights = [
                hypergeom.pmf(np.arange(target + 1), n, count, target)
                for count, n, target in zip(derived, called, projection)
            ]
            result += np.outer(weights[0], weights[1])
            projected += 1
        usable += 1
    if projected == 0:
        raise ValueError("no biallelic records met the projection")
    metadata = {
        "population_names": names,
        "diploid_samples": [len(group) for group in sample_groups],
        "projection": list(projection),
        "biallelic_records": usable,
        "projected_records": projected,
        "polarized": polarized,
        "ancestral_reference_records": ancestral_reference,
        "ancestral_alternate_records": ancestral_alternate,
    }
    return result, metadata


def write_fit(result: dict, path: str | Path) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(result, indent=2) + "\n")
