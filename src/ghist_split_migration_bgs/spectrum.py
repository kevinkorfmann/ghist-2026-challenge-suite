from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Sequence

import numpy as np

from ghist_inference.validation import (
    validate_spectrum,
    validate_fit_settings,
    optimizer_diagnostics,
    validate_expected_spectrum,
)
import pysam

from ghist_inference.vcf import diploid_alleles, validate_position
from scipy.stats import hypergeom

from ghist_growth_bgs.audit import coding_distance, load_coding_intervals


def isolation_with_migration_spectrum(
    parameters: Sequence[float],
    sample_sizes: Sequence[int],
    *,
    folded: bool = False,
    dt_fac: float = 0.02,
    adapt_dt: bool = False,
):
    """Expected joint spectrum under a split with continuous symmetric migration.

    Parameters are ``nu1, nu2, T, M``. Population sizes are relative to the
    ancestral reference, time is in ``2 * Nref`` units, and
    ``M = 2 * Nref * m``.
    """
    import moments

    nu1, nu2, duration, migration = parameters
    if min(nu1, nu2, duration, migration) <= 0:
        raise ValueError("all model parameters must be positive")
    n1, n2 = (int(value) for value in sample_sizes)
    spectrum = moments.Spectrum(moments.LinearSystem_1D.steady_state_1D(n1 + n2))
    spectrum = moments.Manips.split_1D_to_2D(spectrum, n1, n2)
    matrix = np.asarray([[0.0, migration], [migration, 0.0]], dtype=np.float64)
    spectrum.integrate([nu1, nu2], duration, m=matrix, dt_fac=dt_fac, adapt_dt=adapt_dt)
    return spectrum.fold() if folded else spectrum


def _fit_restart(arguments: tuple) -> dict:
    spectrum, start, lower, upper, restart, maxiter, dt_fac, adapt_dt = arguments
    import moments

    data = moments.Spectrum(np.asarray(spectrum, dtype=np.float64))

    def model(parameters, sample_sizes):
        expected = isolation_with_migration_spectrum(
            parameters,
            sample_sizes,
            dt_fac=dt_fac,
            adapt_dt=adapt_dt,
        )
        try:
            return validate_expected_spectrum(expected, data)
        except ValueError:
            if adapt_dt:
                raise
            return validate_expected_spectrum(
                isolation_with_migration_spectrum(
                    parameters,
                    sample_sizes,
                    dt_fac=dt_fac,
                    adapt_dt=True,
                ),
                data,
            )

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


def fit_isolation_with_migration_spectrum(
    spectrum: np.ndarray,
    *,
    mutation_rate: float,
    sequence_length: float,
    projection: tuple[int, int] | None = None,
    restarts: int = 48,
    workers: int = 1,
    maxiter: int = 750,
    seed: int = 75193298,
    integration_dt_fac: float = 0.02,
    adaptive_integration: bool = False,
    initial_scaled: Sequence[float] | None = None,
) -> dict:
    """Fit the four-parameter split-and-migration joint spectrum."""
    import moments

    spectrum = validate_spectrum(spectrum, 2)
    validate_fit_settings(
        sequence_length=sequence_length,
        mutation_rate=mutation_rate,
        restarts=restarts,
        workers=workers,
        maxiter=maxiter,
        shape=spectrum.shape,
        projection=projection,
        integration_dt_fac=integration_dt_fac,
    )
    data = moments.Spectrum(np.asarray(spectrum, dtype=np.float64))
    if projection is not None:
        data = data.project(list(projection))
    lower = np.asarray([0.01, 0.01, 1e-5, 1e-7], dtype=np.float64)
    upper = np.asarray([100.0, 100.0, 20.0, 50.0], dtype=np.float64)
    rng = np.random.default_rng(seed)
    if initial_scaled is None:
        starts = [
            np.exp(rng.uniform(np.log(lower), np.log(upper))) for _ in range(restarts)
        ]
    else:
        initial = np.asarray(initial_scaled, dtype=np.float64)
        if (
            initial.shape != lower.shape
            or not np.isfinite(initial).all()
            or np.any(initial <= lower)
            or np.any(initial >= upper)
        ):
            raise ValueError("initial_scaled must lie strictly within parameter bounds")
        starts = [initial]
        starts.extend(
            np.clip(
                np.exp(np.log(initial) + rng.normal(0.0, 0.5, len(initial))),
                lower * 1.001,
                upper / 1.001,
            )
            for _ in range(restarts - 1)
        )
    work = [
        (
            np.asarray(data),
            start,
            lower,
            upper,
            restart,
            maxiter,
            integration_dt_fac,
            adaptive_integration,
        )
        for restart, start in enumerate(starts)
    ]
    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            results = list(executor.map(_fit_restart, work))
    else:
        results = [_fit_restart(arguments) for arguments in work]
    valid = [row for row in results if row["ok"]]
    if not valid:
        raise RuntimeError(f"all moments restarts failed: {results[:3]}")
    best = max(valid, key=lambda row: row["likelihood"])
    nu1, nu2, duration, scaled_migration = best["estimate"]
    ancestral_size = best["theta"] / (4 * mutation_rate * sequence_length)
    return {
        "AFR_population_size": float(nu1 * ancestral_size),
        "EUR_population_size": float(nu2 * ancestral_size),
        "generations_since_split": float(duration * 2 * ancestral_size),
        "migration_rate": float(scaled_migration / (2 * ancestral_size)),
        "ancestral_population_size": float(ancestral_size),
        "scaled_parameters": {
            "AFR_size_ratio": float(nu1),
            "EUR_size_ratio": float(nu2),
            "split_duration_2N": float(duration),
            "migration_2N": float(scaled_migration),
        },
        "theta": best["theta"],
        "multinomial_log_likelihood": best["likelihood"],
        "optimizer_diagnostics": optimizer_diagnostics(
            best["estimate"],
            lower,
            upper,
            [r["likelihood"] for r in valid],
            best["likelihood"],
        ),
        "best_restart": int(best["restart"]),
        "successful_restarts": len(valid),
        "failed_restarts": [row for row in results if not row["ok"]],
        "projection": [int(value) for value in data.sample_sizes],
        "restarts": restarts,
        "workers": workers,
        "sequence_length": sequence_length,
        "mutation_rate": mutation_rate,
        "integration_dt_fac": integration_dt_fac,
        "adaptive_integration": adaptive_integration,
        "initial_scaled": None if initial_scaled is None else list(initial_scaled),
    }


def stratified_joint_spectra_from_vcf(
    vcf_path: Path,
    coding_bed_path: Path,
) -> tuple[dict[str, np.ndarray], dict]:
    """Extract unfolded AFR/EUR spectra by distance from coding sequence."""
    coding = load_coding_intervals(coding_bed_path)
    starts = {
        contig: [interval[0] for interval in intervals]
        for contig, intervals in coding.items()
    }
    vcf = pysam.VariantFile(str(vcf_path))
    afr = [sample for sample in vcf.header.samples if sample.startswith("AFR_")]
    eur = [sample for sample in vcf.header.samples if sample.startswith("EUR_")]
    if not afr or not eur or len(afr) + len(eur) != len(vcf.header.samples):
        raise ValueError(f"unexpected sample counts: AFR={len(afr)}, EUR={len(eur)}")
    maximum = (2 * len(afr), 2 * len(eur))
    labels = (
        "coding_overlap",
        "noncoding_1bp_1kb",
        "noncoding_1kb_10kb",
        "noncoding_10kb_100kb",
        "noncoding_gt100kb",
    )
    spectra = {
        label: np.zeros((maximum[0] + 1, maximum[1] + 1), dtype=np.int64)
        for label in labels
    }
    counts = {label: 0 for label in labels}
    ancestral_reference = ancestral_alternate = skipped = 0
    lengths = {c: vcf.header.contigs[c].length for c in vcf.header.contigs}
    if set(lengths) != set(coding) or any(
        length is None or length <= 0 for length in lengths.values()
    ):
        raise ValueError(
            "VCF and coding BED need identical contigs with positive lengths"
        )
    for contig, intervals in coding.items():
        if intervals and intervals[-1][1] > lengths[contig]:
            raise ValueError("Coding interval exceeds chromosome length")
    for record in vcf:
        validate_position(record, lengths)
        if len(record.alleles or ()) != 2 or any(len(a) != 1 for a in record.alleles):
            skipped += 1
            continue
        group_counts = []
        group_called = []
        for samples in (afr, eur):
            alleles = diploid_alleles(record, samples, require_complete=True)
            group_called.append(len(alleles))
            group_counts.append(sum(alleles))
        if tuple(group_called) != maximum:
            skipped += 1
            continue
        ancestral = record.info.get("AA")
        if isinstance(ancestral, tuple):
            ancestral = ancestral[0] if ancestral else None
        if ancestral == record.ref:
            ancestral_reference += 1
        elif ancestral == record.alts[0]:
            ancestral_alternate += 1
            group_counts = [n - k for n, k in zip(group_called, group_counts)]
        else:
            skipped += 1
            continue
        distance = coding_distance(
            coding[record.contig], record.pos - 1, starts[record.contig]
        )
        if distance == 0:
            label = "coding_overlap"
        elif distance <= 1_000:
            label = "noncoding_1bp_1kb"
        elif distance <= 10_000:
            label = "noncoding_1kb_10kb"
        elif distance <= 100_000:
            label = "noncoding_10kb_100kb"
        else:
            label = "noncoding_gt100kb"
        spectra[label][tuple(group_counts)] += 1
        counts[label] += 1
    vcf.close()
    metadata = {
        "diploid_samples": {"AFR": len(afr), "EUR": len(eur)},
        "haplotypes": {"AFR": maximum[0], "EUR": maximum[1]},
        "ancestral_reference_records": ancestral_reference,
        "ancestral_alternate_records": ancestral_alternate,
        "skipped_records": skipped,
        "sites_by_distance_stratum": counts,
    }
    return spectra, metadata


def project_joint_spectrum(
    spectrum: np.ndarray, projection: tuple[int, int]
) -> np.ndarray:
    """Project an integer joint spectrum without requiring moments."""
    n1, n2 = np.asarray(spectrum).shape[0] - 1, np.asarray(spectrum).shape[1] - 1
    if not (2 <= projection[0] <= n1 and 2 <= projection[1] <= n2):
        raise ValueError("projection exceeds source dimensions")
    result = np.zeros((projection[0] + 1, projection[1] + 1), dtype=np.float64)
    for i, j in zip(*np.nonzero(spectrum)):
        left = hypergeom.pmf(np.arange(projection[0] + 1), n1, i, projection[0])
        right = hypergeom.pmf(np.arange(projection[1] + 1), n2, j, projection[1])
        result += float(spectrum[i, j]) * np.outer(left, right)
    return result
