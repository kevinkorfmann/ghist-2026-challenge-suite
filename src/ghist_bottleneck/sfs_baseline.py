from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pysam

from ghist_inference.validation import (
    validate_spectrum,
    validate_fit_settings,
    optimizer_diagnostics,
)
from ghist_inference.vcf import diploid_alleles, validate_position
from scipy.stats import hypergeom


def projected_sfs_from_vcf(
    path: str | Path,
    *,
    projection: int | None = None,
    polarized: bool = False,
) -> tuple[np.ndarray, dict]:
    """Build an allele-count spectrum with exact missing-data projection."""
    vcf = pysam.VariantFile(str(path))
    counts: list[tuple[int, int]] = []
    max_called = 0
    ancestral_reference = ancestral_alternate = 0
    lengths = {c: vcf.header.contigs[c].length for c in vcf.header.contigs}
    for rec in vcf:
        validate_position(rec, lengths)
        if len(rec.alleles or ()) != 2 or any(len(a) != 1 for a in rec.alleles):
            continue
        alleles = diploid_alleles(rec, vcf.header.samples)
        called, alt = len(alleles), sum(alleles)
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
            counts.append((derived, called))
            max_called = max(max_called, called)
    vcf.close()
    if not counts:
        raise ValueError("No segregating biallelic SNPs found")
    called_array = np.asarray([n for _, n in counts])
    if projection is None:
        # Keeps at least 90% of variable sites while avoiding a single poor site
        # collapsing the projection size.
        projection = int(np.quantile(called_array, 0.1, method="lower"))
    if projection < 4 or projection > max_called:
        raise ValueError(
            f"invalid projection {projection}; maximum called is {max_called}"
        )
    spectrum = np.zeros(projection + 1, dtype=np.float64)
    kept = 0
    target_counts = np.arange(projection + 1)
    for alt, called in counts:
        if called < projection:
            continue
        spectrum += hypergeom.pmf(target_counts, called, alt, projection)
        kept += 1
    metadata = {
        "projection": projection,
        "segregating_records": len(counts),
        "projected_records": kept,
        "maximum_called_haplotypes": max_called,
        "polarized": polarized,
        "ancestral_reference_records": ancestral_reference,
        "ancestral_alternate_records": ancestral_alternate,
    }
    return spectrum, metadata


def fit_two_epoch(
    path: str | Path,
    *,
    mutation_rate: float = 3.42e-8,
    sequence_length: float = 100_000_000,
    restarts: int = 24,
    seed: int = 4721,
    polarized: bool = True,
    folded: bool = False,
) -> dict:
    """Fit the mechanistic two-epoch SFS reference with moments."""
    spectrum, metadata = projected_sfs_from_vcf(path, polarized=polarized)
    result = fit_two_epoch_spectrum(
        spectrum,
        mutation_rate=mutation_rate,
        sequence_length=sequence_length,
        restarts=restarts,
        seed=seed,
        folded=folded,
    )
    return {**result, "polarized": polarized, "folded": folded, "sfs": metadata}


def fit_two_epoch_spectrum(
    spectrum: np.ndarray,
    *,
    mutation_rate: float = 3.42e-8,
    sequence_length: float = 100_000_000,
    restarts: int = 24,
    seed: int = 4721,
    folded: bool = False,
    initial: tuple[float, float] | None = None,
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
        return expectation.fold() if folded else expectation

    lower = [0.02, 1e-4]
    upper = [0.999, 5.0]
    rng = np.random.default_rng(seed)
    fits = []
    failures = []
    if initial is not None:
        initial = np.asarray(initial, dtype=float)
        if (
            initial.shape != (2,)
            or not np.isfinite(initial).all()
            or np.any(initial <= lower)
            or np.any(initial >= upper)
        ):
            raise ValueError(
                "Initial estimate must lie strictly within parameter bounds"
            )
    for restart in range(restarts):
        if restart == 0 and initial is not None:
            start = np.asarray(initial, dtype=float)
        else:
            start = np.exp(rng.uniform(np.log(lower), np.log(upper)))
        try:
            estimate = moments.Inference.optimize_log_lbfgsb(
                start,
                data,
                model,
                lower_bound=lower,
                upper_bound=upper,
                maxiter=500,
                verbose=0,
            )
            expected = model(estimate, data.sample_sizes)
            likelihood = float(moments.Inference.ll_multinom(expected, data))
            theta = float(moments.Inference.optimal_sfs_scaling(expected, data))
            if not np.isfinite([*estimate, likelihood, theta]).all() or theta <= 0:
                raise ValueError("Nonfinite fit or nonpositive theta")
            fits.append((likelihood, np.asarray(estimate), theta))
        except Exception as error:
            failures.append(
                {"restart": restart, "reason": f"{type(error).__name__}: {error}"}
            )
    if not fits:
        raise RuntimeError(f"All moments restarts failed: {failures[:3]}")
    likelihood, estimate, theta = max(fits, key=lambda item: item[0])
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
        "successful_restarts": len(fits),
        "failed_restarts": failures,
        "optimizer_diagnostics": optimizer_diagnostics(
            estimate, lower, upper, [fit[0] for fit in fits], likelihood
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
