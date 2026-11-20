from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pysam


@dataclass(frozen=True)
class VariantMatrix:
    contig: str
    contig_length: int
    positions: np.ndarray
    haplotypes: np.ndarray
    skipped_records: int


def read_phased_biallelic_vcf(path: Path) -> VariantMatrix:
    positions: list[int] = []
    rows: list[list[int]] = []
    skipped = 0
    phase_sets: dict[str, int | None] = {}
    with pysam.VariantFile(str(path)) as vcf:
        contigs = list(vcf.header.contigs)
        if len(contigs) != 1:
            raise ValueError("expected exactly one contig")
        contig = contigs[0]
        contig_length = int(vcf.header.contigs[contig].length)
        previous_position = -1
        for record in vcf:
            position = int(record.pos) - 1
            if record.contig != contig or not 0 <= position < contig_length:
                raise ValueError(
                    "VCF coordinates must be 1-based and within the declared contig"
                )
            if position <= previous_position:
                raise ValueError("VCF SNP positions must be sorted and unique")
            previous_position = position
            if len(record.alleles or ()) != 2:
                skipped += 1
                continue
            if any(
                len(allele) != 1 or allele.upper() not in "ACGT"
                for allele in record.alleles
            ):
                raise ValueError(
                    "Sweep inference requires single-nucleotide A/C/G/T alleles"
                )
            row: list[int] = []
            valid = True
            row_phase_sets = {}
            for sample_name, sample in record.samples.items():
                gt = sample.get("GT")
                if gt is None or len(gt) != 2 or None in gt:
                    raise ValueError(
                        "Missing or non-diploid genotypes require a callable-exposure model"
                    )
                if not sample.phased and gt[0] != gt[1]:
                    raise ValueError(
                        "Haplotype sweep inference requires phased heterozygous genotypes"
                    )
                if gt[0] > 1 or gt[1] > 1:
                    valid = False
                    break
                if gt[0] != gt[1]:
                    phase_set = sample.get("PS")
                    if phase_set is not None and (
                        type(phase_set) is not int or not 0 <= phase_set <= 2**31 - 1
                    ):
                        raise ValueError(
                            "PS must be a nonnegative 32-bit integer or missing"
                        )
                    row_phase_sets[sample_name] = phase_set
                row.extend((int(gt[0]), int(gt[1])))
            if not valid:
                skipped += 1
                continue
            if sum(row) in (0, len(row)):
                skipped += 1
                continue
            # Homozygotes and skipped records cannot connect phase blocks.
            # Missing PS denotes one implicit chromosome-wide phase set.
            for sample_name, phase_set in row_phase_sets.items():
                if sample_name in phase_sets and phase_sets[sample_name] != phase_set:
                    raise ValueError(
                        f"Sample {sample_name} has mixed or multiple PS phase sets; sweep LD requires consistent chromosome-wide phase"
                    )
                phase_sets[sample_name] = phase_set
            positions.append(position)
            rows.append(row)
    matrix = np.asarray(rows, dtype=np.uint8)
    if matrix.ndim != 2 or matrix.shape[1] < 4:
        raise ValueError("VCF did not yield a usable haplotype matrix")
    return VariantMatrix(
        contig=contig,
        contig_length=contig_length,
        positions=np.asarray(positions, dtype=np.int32),
        haplotypes=matrix,
        skipped_records=skipped,
    )


def read_bed_intervals(path: Path, contig: str) -> np.ndarray:
    intervals: list[tuple[int, int]] = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        fields = line.split()
        if fields[0].removeprefix("chr") == contig.removeprefix("chr"):
            intervals.append((int(fields[1]), int(fields[2])))
    if not intervals:
        raise ValueError(f"no accessibility intervals for {contig}")
    return merge_intervals(np.asarray(intervals, dtype=np.int64))


def merge_intervals(intervals: np.ndarray) -> np.ndarray:
    intervals = np.asarray(intervals, dtype=np.int64)
    if intervals.ndim != 2 or intervals.shape[1] != 2:
        raise ValueError("Expected an array of BED interval pairs")
    merged = []
    for left, right in sorted(intervals.tolist()):
        if left < 0 or right <= left:
            raise ValueError("Invalid half-open accessibility interval")
        if merged and left <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], right)
        else:
            merged.append([left, right])
    return np.asarray(merged, dtype=np.int64).reshape(-1, 2)


def accessible_bases(
    starts: np.ndarray, ends: np.ndarray, intervals: np.ndarray
) -> np.ndarray:
    result = np.zeros(starts.size, dtype=np.int64)
    for left, right in merge_intervals(intervals):
        result += np.maximum(0, np.minimum(ends, right) - np.maximum(starts, left))
    return result


def tajima_d(
    allele_counts: np.ndarray, sequence_length: int, sample_size: int
) -> float:
    allele_counts = np.asarray(allele_counts)
    allele_counts = allele_counts[(allele_counts > 0) & (allele_counts < sample_size)]
    segregating = int(allele_counts.size)
    if segregating < 2 or sequence_length <= 0:
        return float("nan")
    a1 = np.sum(1.0 / np.arange(1, sample_size))
    a2 = np.sum(1.0 / np.arange(1, sample_size) ** 2)
    b1 = (sample_size + 1.0) / (3.0 * (sample_size - 1.0))
    b2 = (
        2.0
        * (sample_size**2 + sample_size + 3.0)
        / (9.0 * sample_size * (sample_size - 1.0))
    )
    c1 = b1 - 1.0 / a1
    c2 = b2 - (sample_size + 2.0) / (a1 * sample_size) + a2 / a1**2
    e1 = c1 / a1
    e2 = c2 / (a1**2 + a2)
    pi_total = np.sum(2.0 * allele_counts * (sample_size - allele_counts)) / (
        sample_size * (sample_size - 1.0)
    )
    denominator = np.sqrt(e1 * segregating + e2 * segregating * (segregating - 1.0))
    return (
        float((pi_total - segregating / a1) / denominator)
        if denominator > 0
        else float("nan")
    )


def haplotype_homozygosity(window_haplotypes: np.ndarray) -> tuple[float, float, float]:
    if window_haplotypes.shape[0] == 0:
        return float("nan"), float("nan"), float("nan")
    _, counts = np.unique(window_haplotypes.T, axis=0, return_counts=True)
    frequencies = np.sort(counts.astype(float) / window_haplotypes.shape[1])[::-1]
    h1 = float(np.sum(frequencies**2))
    h12 = float((frequencies[:2].sum()) ** 2 + np.sum(frequencies[2:] ** 2))
    h2_h1 = float((h1 - frequencies[0] ** 2) / h1) if h1 > 0 else float("nan")
    return h1, h12, h2_h1


def mean_pairwise_r2(window_haplotypes: np.ndarray, maximum_sites: int = 96) -> float:
    if window_haplotypes.shape[0] < 2:
        return float("nan")
    allele_counts = window_haplotypes.sum(axis=1)
    keep = (allele_counts >= 2) & (allele_counts <= window_haplotypes.shape[1] - 2)
    values = window_haplotypes[keep].astype(float)
    if values.shape[0] < 2:
        return float("nan")
    if values.shape[0] > maximum_sites:
        indices = np.linspace(0, values.shape[0] - 1, maximum_sites).round().astype(int)
        values = values[indices]
    values -= values.mean(axis=1, keepdims=True)
    scales = np.sqrt(np.sum(values**2, axis=1, keepdims=True))
    values = values / np.maximum(scales, 1e-12)
    correlation = values @ values.T
    upper = correlation[np.triu_indices(correlation.shape[0], 1)]
    return float(np.mean(upper**2))


def robust_z(values: np.ndarray) -> np.ndarray:
    result = np.full(values.shape, np.nan, dtype=float)
    valid = np.isfinite(values)
    if not np.any(valid):
        return result
    median = np.median(values[valid])
    mad = np.median(np.abs(values[valid] - median))
    scale = max(1.4826 * mad, 1e-12)
    result[valid] = (values[valid] - median) / scale
    return result


def scan_windows(
    variants: VariantMatrix,
    accessibility: np.ndarray,
    window_size: int,
    step_size: int,
) -> dict[str, np.ndarray]:
    if window_size <= 0 or step_size <= 0:
        raise ValueError("Window and step sizes must be positive")
    intervals = merge_intervals(accessibility)
    accessible = np.zeros(variants.positions.size, dtype=bool)
    if len(intervals):
        # Avoid comparing every variant with every interval on chromosome-scale data.
        interval_index = (
            np.searchsorted(intervals[:, 0], variants.positions, side="right") - 1
        )
        accessible = (interval_index >= 0) & (
            variants.positions < intervals[np.maximum(interval_index, 0), 1]
        )
    positions = variants.positions[accessible]
    variant_haplotypes = variants.haplotypes[accessible]
    starts = np.arange(0, variants.contig_length, step_size, dtype=np.int64)
    ends = np.minimum(starts + window_size, variants.contig_length)
    exposure = accessible_bases(starts, ends, accessibility)
    n_windows = starts.size
    features = {
        "site_count": np.zeros(n_windows, dtype=float),
        "pi": np.full(n_windows, np.nan),
        "tajima_d": np.full(n_windows, np.nan),
        "singleton_fraction": np.full(n_windows, np.nan),
        "intermediate_fraction": np.full(n_windows, np.nan),
        "h1": np.full(n_windows, np.nan),
        "h12": np.full(n_windows, np.nan),
        "h2_h1": np.full(n_windows, np.nan),
        "mean_r2": np.full(n_windows, np.nan),
    }
    sample_size = variants.haplotypes.shape[1]
    for index, (left, right) in enumerate(zip(starts, ends, strict=True)):
        lo = int(np.searchsorted(positions, left, side="left"))
        hi = int(np.searchsorted(positions, right, side="left"))
        haplotypes = variant_haplotypes[lo:hi]
        if haplotypes.shape[0] == 0 or exposure[index] == 0:
            continue
        counts = haplotypes.sum(axis=1).astype(float)
        minor = np.minimum(counts, sample_size - counts)
        features["site_count"][index] = haplotypes.shape[0]
        pi_total = np.sum(2.0 * counts * (sample_size - counts)) / (
            sample_size * (sample_size - 1.0)
        )
        features["pi"][index] = pi_total / exposure[index]
        features["tajima_d"][index] = tajima_d(
            counts, int(exposure[index]), sample_size
        )
        features["singleton_fraction"][index] = np.mean(minor == 1)
        features["intermediate_fraction"][index] = np.mean(minor >= sample_size * 0.2)
        h1, h12, h2_h1 = haplotype_homozygosity(haplotypes)
        features["h1"][index] = h1
        features["h12"][index] = h12
        features["h2_h1"][index] = h2_h1
        features["mean_r2"][index] = mean_pairwise_r2(haplotypes)

    span = ends - starts
    valid_exposure = (span >= window_size * 0.8) & (exposure >= span * 0.8)
    # Uncallable and truncated windows must not set the chromosome's reference
    # distribution. Masking only the final score still changes all other scores.
    normalized = {
        key: np.where(valid_exposure, values, np.nan)
        for key, values in features.items()
    }
    components = {
        "low_pi": -robust_z(np.log(np.maximum(normalized["pi"], 1e-12))),
        "low_sites": -robust_z(
            np.log1p(normalized["site_count"] / np.maximum(exposure, 1))
        ),
        "low_tajima_d": -robust_z(normalized["tajima_d"]),
        "high_h12": robust_z(np.log(np.maximum(normalized["h12"], 1e-6))),
        "high_h2_h1": robust_z(normalized["h2_h1"]),
        "high_r2": robust_z(np.log(np.maximum(normalized["mean_r2"], 1e-6))),
    }
    components["conditioned_h2_h1"] = np.where(
        components["high_h12"] > 1.0,
        np.maximum(components["high_h2_h1"], 0.0),
        0.0,
    )
    hard_component_names = [
        "low_pi",
        "low_sites",
        "low_tajima_d",
        "high_h12",
        "high_r2",
    ]
    stack = np.vstack(
        [np.clip(components[name], -5.0, 8.0) for name in hard_component_names]
    )
    counts = np.isfinite(stack).sum(axis=0)
    score = np.divide(
        np.nansum(stack, axis=0),
        counts,
        out=np.full(starts.shape, np.nan),
        where=counts > 0,
    )
    score[~valid_exposure] = np.nan
    soft_weights = {
        "low_pi": 0.5,
        "low_tajima_d": 0.75,
        "high_h12": 1.5,
        "conditioned_h2_h1": 0.5,
        "high_r2": 0.5,
    }
    soft_stack = np.vstack(
        [
            weight * np.clip(components[name], -5.0, 8.0)
            for name, weight in soft_weights.items()
        ]
    )
    soft_score = np.nansum(soft_stack, axis=0) / sum(soft_weights.values())
    soft_score[~valid_exposure] = np.nan
    return {
        "start": starts,
        "end": ends,
        "center": (starts + ends) // 2,
        "accessible_bases": exposure,
        **features,
        **{f"z_{key}": value for key, value in components.items()},
        "composite_score": score,
        "soft_composite_score": soft_score,
    }


def flank_ld_scan(
    variants: VariantMatrix,
    candidate_positions: np.ndarray,
    flank_width: int,
    maximum_sites_per_flank: int = 128,
) -> dict[str, np.ndarray]:
    sample_size = variants.haplotypes.shape[1]
    scores = np.full(candidate_positions.size, np.nan)
    within_left = np.full(candidate_positions.size, np.nan)
    within_right = np.full(candidate_positions.size, np.nan)
    between = np.full(candidate_positions.size, np.nan)

    for index, center in enumerate(candidate_positions):
        bounds = [center - flank_width, center, center + flank_width]
        lo = int(np.searchsorted(variants.positions, bounds[0], side="left"))
        mid = int(np.searchsorted(variants.positions, bounds[1], side="left"))
        hi = int(np.searchsorted(variants.positions, bounds[2], side="left"))
        left = variants.haplotypes[lo:mid]
        right = variants.haplotypes[mid:hi]
        prepared: list[np.ndarray] = []
        for values in (left, right):
            allele_counts = values.sum(axis=1)
            keep = (allele_counts >= 2) & (allele_counts <= sample_size - 2)
            values = values[keep].astype(float)
            if values.shape[0] > maximum_sites_per_flank:
                selected = np.linspace(0, values.shape[0] - 1, maximum_sites_per_flank)
                values = values[np.rint(selected).astype(int)]
            if values.shape[0] < 4:
                prepared = []
                break
            values -= values.mean(axis=1, keepdims=True)
            values /= np.maximum(
                np.sqrt(np.sum(values**2, axis=1, keepdims=True)), 1e-12
            )
            prepared.append(values)
        if len(prepared) != 2:
            continue
        left_std, right_std = prepared
        left_corr = left_std @ left_std.T
        right_corr = right_std @ right_std.T
        cross_corr = left_std @ right_std.T
        left_r2 = float(np.mean(left_corr[np.triu_indices(left_corr.shape[0], 1)] ** 2))
        right_r2 = float(
            np.mean(right_corr[np.triu_indices(right_corr.shape[0], 1)] ** 2)
        )
        cross_r2 = float(np.mean(cross_corr**2))
        within_left[index] = left_r2
        within_right[index] = right_r2
        between[index] = cross_r2
        scores[index] = np.sqrt(left_r2 * right_r2) / max(cross_r2, 1e-12)
    return {
        "candidate_position": candidate_positions,
        "flank_width": np.full(candidate_positions.size, flank_width),
        "within_left_r2": within_left,
        "within_right_r2": within_right,
        "between_r2": between,
        "flank_ld_ratio": scores,
    }
