from __future__ import annotations

import json
import hashlib
from collections import Counter
from pathlib import Path

import numpy as np
import pysam

from ghist_inference.validation import (
    validate_spectrum,
    validate_fit_settings,
    optimizer_diagnostics,
    validate_expected_spectrum,
)
from ghist_inference.vcf import diploid_alleles, validate_position
from ghist_inference.checkpoint import spectrum_digest
from scipy.stats import hypergeom
from .optimization import DEFAULT_STRATEGY, LEGACY_STRATEGY, LOWER, UPPER, log_starts, run_optimizer


def projected_sfs_from_vcf(
    path: str | Path,
    *,
    projection: int | None = None,
    polarized: bool = False,
) -> tuple[np.ndarray, dict]:
    """Build an allele-count spectrum with exact missing-data projection."""
    vcf = pysam.VariantFile(str(path))
    counts: Counter[tuple[int, int]] = Counter()
    max_called = 0
    ancestral_reference = ancestral_alternate = 0
    incomplete_records = 0
    lengths = {c: vcf.header.contigs[c].length for c in vcf.header.contigs}
    for rec in vcf:
        validate_position(rec, lengths)
        if len(rec.alleles or ()) != 2 or any(len(a) != 1 for a in rec.alleles):
            continue
        alleles = diploid_alleles(rec, vcf.header.samples)
        called, alt = len(alleles), sum(alleles)
        incomplete_records += called < 2 * len(vcf.header.samples)
        derived = alt
        if polarized:
            ancestral = rec.info.get("AA")
            if isinstance(ancestral, tuple):
                ancestral = ancestral[0] if ancestral else None
            if ancestral == rec.ref:
                ancestral_reference += 1
            elif ancestral == rec.alts[0]:
                ancestral_alternate += 1
                derived = called - alt
            else:
                raise ValueError(
                    f"record {rec.contig}:{rec.pos} has invalid ancestral allele {ancestral!r}"
                )
        if 0 < derived < called:
            counts[derived, called] += 1
            max_called = max(max_called, called)
    vcf.close()
    if not counts:
        raise ValueError("No segregating biallelic SNPs found")
    total_records = sum(counts.values())
    if projection is None:
        # Keeps at least 90% of variable sites while avoiding a single poor site
        # collapsing the projection size.
        by_called = Counter()
        for (_, called), multiplicity in counts.items():
            by_called[called] += multiplicity
        rank = int(np.floor(.1 * (total_records - 1)))
        cumulative = 0
        for called, multiplicity in sorted(by_called.items()):
            cumulative += multiplicity
            if cumulative > rank:
                projection = called
                break
    if projection < 4 or projection > max_called:
        raise ValueError(
            f"invalid projection {projection}; maximum called is {max_called}"
        )
    spectrum = np.zeros(projection + 1, dtype=np.float64)
    kept = 0
    target_counts = np.arange(projection + 1)
    for (alt, called), multiplicity in counts.items():
        if called < projection:
            continue
        if called == projection:
            spectrum[alt] += multiplicity
        else:
            spectrum += multiplicity * hypergeom.pmf(target_counts, called, alt, projection)
        kept += multiplicity
    metadata = {
        "spectrum_sha256": spectrum_digest(spectrum),
        "projection": projection,
        "segregating_records": total_records,
        "projected_records": kept,
        "maximum_called_haplotypes": max_called,
        "incomplete_genotype_records": incomplete_records,
        "polarized": polarized,
        "ancestral_reference_records": ancestral_reference,
        "ancestral_alternate_records": ancestral_alternate,
    }
    return spectrum, metadata


def spectrum_exposure(metadata: dict, sequence_length: float,
                      callable_sequence_length: float | None) -> float:
    """Do not infer callable exposure from variant-only genotype calls."""
    if callable_sequence_length is None:
        if metadata.get("incomplete_genotype_records", 0):
            raise ValueError(
                "Incomplete genotypes require an independently measured "
                "callable_sequence_length for the chosen projection; a "
                "variant-only VCF cannot establish callable exposure."
            )
        return sequence_length
    if not np.isfinite(callable_sequence_length) or not 0 < callable_sequence_length <= sequence_length:
        raise ValueError("callable_sequence_length must be positive and within sequence_length")
    return callable_sequence_length


def fit_two_epoch(
    path: str | Path,
    *,
    mutation_rate: float = 3.42e-8,
    sequence_length: float = 100_000_000,
    callable_sequence_length: float | None = None,
    restarts: int = 24,
    seed: int = 4721,
    polarized: bool = True,
    folded: bool = False,
    optimizer_strategy: str = DEFAULT_STRATEGY,
) -> dict:
    """Fit the mechanistic two-epoch SFS reference with moments."""
    if not polarized and not folded:
        raise ValueError('Unpolarized VCF input requires a folded spectrum')
    with Path(path).open("rb") as handle:
        input_sha256 = hashlib.file_digest(handle, "sha256").hexdigest()
    spectrum, metadata = projected_sfs_from_vcf(path, polarized=polarized)
    exposure = spectrum_exposure(metadata, sequence_length, callable_sequence_length)
    result = fit_two_epoch_spectrum(
        spectrum,
        mutation_rate=mutation_rate,
        sequence_length=exposure,
        restarts=restarts,
        seed=seed,
        folded=folded,
        optimizer_strategy=optimizer_strategy,
    )
    with Path(path).open('rb') as handle:
        if hashlib.file_digest(handle, 'sha256').hexdigest() != input_sha256:
            raise ValueError("VCF changed during spectrum extraction or fitting; regenerate the fit")
    return {**result, "polarized": polarized, "folded": folded, "sfs": metadata,
            "input_vcf_sha256": input_sha256,
            "sequence_length": sequence_length, "callable_sequence_length": exposure}


def fit_two_epoch_spectrum(
    spectrum: np.ndarray,
    *,
    mutation_rate: float = 3.42e-8,
    sequence_length: float = 100_000_000,
    restarts: int = 24,
    seed: int = 4721,
    folded: bool = False,
    initial: tuple[float, float] | None = None,
    optimizer_strategy: str = DEFAULT_STRATEGY,
    initial_log_starts: np.ndarray | None = None,
) -> dict:
    """Fit a two-epoch model to an already projected allele-count spectrum."""
    import moments

    spectrum = validate_spectrum(spectrum, 1)
    validate_fit_settings(
        sequence_length=sequence_length,
        mutation_rate=mutation_rate,
        restarts=restarts,
        workers=1,
        maxiter=500,
        shape=spectrum.shape,
    )
    data = moments.Spectrum(spectrum)
    if folded:
        data = data.fold()

    def model(params, ns):
        expectation = moments.Demographics1D.two_epoch(params, ns)
        expectation = expectation.fold() if folded else expectation
        return validate_expected_spectrum(expectation, data)

    def objective(point):
        value = -float(moments.Inference.ll_multinom(model(np.exp(point), data.sample_sizes), data))
        if not np.isfinite(value):
            raise ValueError("Nonfinite fixed-domain objective")
        return value

    def evaluate(parameters):
        theta = float(moments.Inference.optimal_sfs_scaling(model(parameters, data.sample_sizes), data))
        if not np.isfinite(theta) or theta <= 0:
            raise ValueError("Nonfinite or nonpositive theta")
        return {"theta": theta}

    def legacy_optimize(point):
        estimate, _, optimizer = moments.Inference.optimize_log_lbfgsb(
            np.exp(point), data, model, lower_bound=LOWER.tolist(), upper_bound=UPPER.tolist(),
            maxiter=500, verbose=0, full_output=True,
        )
        return np.log(estimate), optimizer

    starts = log_starts(restarts, seed, initial=initial, starts=initial_log_starts)
    optimization = run_optimizer(objective, engine="moments", starts=starts, strategy=optimizer_strategy,
                                 evaluate=evaluate, legacy_optimize=legacy_optimize if optimizer_strategy == LEGACY_STRATEGY else None)
    best = optimization["best"]
    estimate = np.asarray(best["parameters"])
    likelihood, theta, best_restart = -best["objective"], best["theta"], best["restart"]
    optimizer = best["optimizer"]
    fraction, scaled_time = estimate
    ancestral_ne = theta / (4 * mutation_rate * sequence_length)
    generations = scaled_time * 2 * ancestral_ne
    return {
        "generations": float(generations),
        "post_decline_fraction": float(fraction),
        "ancestral_effective_size": float(ancestral_ne),
        "scaled_time_2N": float(scaled_time),
        "theta": theta,
        "multinomial_log_likelihood": likelihood,
        "restarts": restarts,
        "mutation_rate": mutation_rate,
        "sequence_length": sequence_length,
        "folded": folded,
        "best_restart": best_restart,
        "engine": "moments",
        "optimizer_contract": optimization["contract"],
        "sample_haplotypes": int(data.sample_sizes[0]),
        "optimizer_strategy": optimizer_strategy,
        "best_optimizer_arm": best["arm"],
        "best_restart_converged": optimizer["converged"],
        "converged_restarts": sum(row["optimizer"]["converged"] for row in optimization["origin_best"]),
        "optimizer": optimizer,
        "optimizer_neighborhood": best["neighborhood"],
        "optimizer_attempts": optimization["attempts"],
        "optimizer_arm_summary": optimization["arm_summary"],
        "optimizer_seed": seed,
        "optimizer_initialization": "explicit_log_starts" if initial_log_starts is not None else "seeded_random_with_optional_first_initial",
        "optimizer_origin_best": [{"restart": row["restart"], "arm": row["arm"], "objective": row["objective"],
                                   "converged": row["optimizer"]["converged"]} for row in optimization["origin_best"]],
        "successful_restarts": len(optimization["origin_best"]),
        "failed_restarts": optimization["failed_origins"],
        "optimizer_diagnostics": optimizer_diagnostics(
            estimate, LOWER, UPPER, [-row["objective"] for row in optimization["origin_best"]], likelihood
        ),
    }


def complete_block_spectra_from_vcf(
    path: str | Path,
    *,
    sequence_length: int = 100_000_000,
    block_length: int = 2_000_000,
    polarized: bool = True,
) -> tuple[np.ndarray, dict]:
    """Extract exact per-block spectra from a complete biallelic challenge VCF."""
    if sequence_length <= 0 or block_length <= 0 or sequence_length % block_length:
        raise ValueError(
            "positive, evenly dividing sequence and block lengths are required"
        )
    vcf = pysam.VariantFile(str(path))
    samples = list(vcf.header.samples)
    projection = 2 * len(samples)
    if projection < 2:
        raise ValueError("VCF requires diploid samples")
    spectra = np.zeros(
        (sequence_length // block_length, projection + 1), dtype=np.float64
    )
    contig = None
    records = np.zeros(sequence_length // block_length, dtype=np.int64)
    for rec in vcf:
        if contig is None:
            contig = rec.contig
        if rec.contig != contig:
            raise ValueError("Block bootstrap expects one contiguous chromosome")
        validate_position(rec, {contig: sequence_length})
        if len(rec.alleles or ()) != 2 or any(
            len(allele) != 1 for allele in rec.alleles
        ):
            continue
        alleles = diploid_alleles(rec, samples, require_complete=True)
        called, alt = len(alleles), sum(alleles)
        derived = alt
        if polarized:
            ancestral = rec.info.get("AA")
            if isinstance(ancestral, tuple):
                ancestral = ancestral[0] if ancestral else None
            if ancestral == rec.alts[0]:
                derived = called - alt
            elif ancestral != rec.ref:
                raise ValueError(f"invalid ancestral allele at {rec.contig}:{rec.pos}")
        if 0 < derived < called:
            block = (rec.pos - 1) // block_length
            if block >= len(records):
                raise ValueError(
                    f"record {rec.contig}:{rec.pos} exceeds the sequence length"
                )
            spectra[block, derived] += 1
            records[block] += 1
    vcf.close()
    if not records.sum():
        raise ValueError("no usable variants found")
    return spectra, {
        "projection": int(projection),
        "blocks": int(len(records)),
        "block_length": int(block_length),
        "records_per_block": records.tolist(),
        "polarized": polarized,
    }


def write_fit(result: dict, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2) + "\n")
