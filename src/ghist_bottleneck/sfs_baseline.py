from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pysam
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
    for rec in vcf:
        if len(rec.alleles or ()) != 2 or any(len(a) != 1 for a in rec.alleles):
            continue
        alt = called = 0
        for sample in vcf.header.samples:
            gt = rec.samples[sample].get("GT")
            if gt is None:
                continue
            vals = [a for a in gt if a is not None and a >= 0]
            called += len(vals)
            alt += sum(a == 1 for a in vals)
        derived = alt
        if polarized:
            ancestral = rec.info.get("AA")
            if isinstance(ancestral, tuple):
                ancestral = ancestral[0]
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
    if not counts:
        raise ValueError("No segregating biallelic SNPs found")
    called_array = np.asarray([n for _, n in counts])
    if projection is None:
        # Keeps at least 90% of variable sites while avoiding a single poor site
        # collapsing the projection size.
        projection = int(np.quantile(called_array, 0.1, method="lower"))
    if projection < 4 or projection > max_called:
        raise ValueError(f"invalid projection {projection}; maximum called is {max_called}")
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
    import moments

    spectrum, metadata = projected_sfs_from_vcf(path, polarized=polarized)
    data = moments.Spectrum(spectrum)
    if folded:
        data = data.fold()

    def model(params, ns):
        expectation = moments.Demographics1D.two_epoch(params, ns)
        return expectation.fold() if folded else expectation

    lower = [0.02, 1e-4]
    upper = [2.0, 5.0]
    rng = np.random.default_rng(seed)
    fits = []
    for _ in range(restarts):
        initial = np.exp(rng.uniform(np.log(lower), np.log(upper)))
        estimate = moments.Inference.optimize_log_lbfgsb(
            initial,
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
        fits.append((likelihood, np.asarray(estimate), theta))
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
        "polarized": polarized,
        "folded": folded,
        "sfs": metadata,
    }


def write_fit(result: dict, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2) + "\n")
