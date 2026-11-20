from __future__ import annotations

import bisect
import collections
import hashlib
from pathlib import Path

import numpy as np
import pysam

from ghist_inference.vcf import diploid_alleles, validate_position


WINDOW_SIZE = 100_000


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    merged: list[list[int]] = []
    for start, end in sorted(intervals):
        if start < 0 or end <= start:
            raise ValueError(f"invalid half-open interval [{start}, {end})")
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(start, end) for start, end in merged]


def load_coding_intervals(path: Path) -> dict[str, list[tuple[int, int]]]:
    raw: dict[str, list[tuple[int, int]]] = collections.defaultdict(list)
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 3:
                raise ValueError(f"BED line {line_number} has fewer than three columns")
            raw[fields[0]].append((int(fields[1]), int(fields[2])))
    return {contig: merge_intervals(intervals) for contig, intervals in raw.items()}


def load_recombination_map(path: Path) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    rows: dict[str, list[tuple[int, float]]] = collections.defaultdict(list)
    with path.open() as handle:
        header = handle.readline().strip().split("\t")
        if header[:3] != ["Chromosome", "Position(bp)", "Rate(cM/Mb)"]:
            raise ValueError(f"unexpected recombination-map header: {header}")
        for line_number, line in enumerate(handle, start=2):
            if not line.strip():
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 3:
                raise ValueError(f"map line {line_number} has fewer than three columns")
            rows[fields[0]].append((int(fields[1]), float(fields[2])))
    result = {}
    for contig, values in rows.items():
        positions = np.asarray([value[0] for value in values], dtype=np.int64)
        rates = np.asarray([value[1] for value in values], dtype=np.float64)
        if np.any(np.diff(positions) <= 0):
            raise ValueError(f"map positions are not strictly increasing on {contig}")
        if np.any(rates < 0) or not np.all(np.isfinite(rates)):
            raise ValueError(f"map rates are invalid on {contig}")
        result[contig] = positions, rates
    return result


def coding_distance(
    intervals: list[tuple[int, int]],
    position: int,
    starts: list[int] | None = None,
) -> int | None:
    if not intervals:
        return None
    if starts is None:
        starts = [interval[0] for interval in intervals]
    index = bisect.bisect_right(starts, position) - 1
    if index >= 0 and position < intervals[index][1]:
        return 0
    distances = []
    if index >= 0:
        distances.append(position - intervals[index][1] + 1)
    if index + 1 < len(intervals):
        distances.append(intervals[index + 1][0] - position)
    return max(0, min(distances)) if distances else None


def local_recombination_rate(
    recombination: tuple[np.ndarray, np.ndarray] | None, position: int
) -> float | None:
    if recombination is None:
        return None
    positions, rates = recombination
    index = int(np.searchsorted(positions, position, side="right") - 1)
    if index < 0:
        return None
    return float(rates[index])


def distance_label(distance: int | None) -> str:
    if distance is None:
        return "unannotated_contig"
    if distance == 0:
        return "coding_overlap"
    if distance <= 1_000:
        return "noncoding_1bp_1kb"
    if distance <= 10_000:
        return "noncoding_1kb_10kb"
    if distance <= 100_000:
        return "noncoding_10kb_100kb"
    return "noncoding_gt100kb"


def overlap_length(intervals: list[tuple[int, int]], start: int, end: int) -> int:
    total = 0
    for left, right in intervals:
        if left >= end:
            break
        if right > start:
            total += max(0, min(right, end) - max(left, start))
    return total


def physical_stratum_lengths(
    contig_lengths: dict[str, int],
    coding: dict[str, list[tuple[int, int]]],
    recombination: dict[str, tuple[np.ndarray, np.ndarray]],
    rate_label,
) -> tuple[dict[str, int], dict[str, int], dict[str, int]]:
    distance_bp: collections.Counter[str] = collections.Counter()
    rate_bp: collections.Counter[str] = collections.Counter()
    joint_bp: collections.Counter[str] = collections.Counter()
    thresholds = (1_000, 10_000, 100_000)
    for contig, length in contig_lengths.items():
        intervals = coding[contig]
        starts = [interval[0] for interval in intervals]
        breakpoints = {0, length}
        for start, end in intervals:
            breakpoints.update((start, end))
            for threshold in thresholds:
                breakpoints.add(max(0, start - threshold))
                breakpoints.add(min(length, end + threshold))
        positions, _ = recombination[contig]
        breakpoints.update(
            int(position) for position in positions if 0 < position < length
        )
        ordered = sorted(breakpoints)
        for left, right in zip(ordered[:-1], ordered[1:]):
            if right <= left:
                continue
            dlabel = distance_label(coding_distance(intervals, left, starts))
            rlabel = rate_label(local_recombination_rate(recombination[contig], left))
            span = right - left
            distance_bp[dlabel] += span
            rate_bp[rlabel] += span
            joint_bp[f"{dlabel}|{rlabel}"] += span
    return dict(distance_bp), dict(rate_bp), dict(joint_bp)


def _empty_stratum() -> dict:
    return {
        "biallelic_sites": 0,
        "oriented_sites": 0,
        "folded_sfs": collections.Counter(),
        "unfolded_sfs": collections.Counter(),
    }


def _add_site(stratum: dict, folded: int, derived: int | None) -> None:
    stratum["biallelic_sites"] += 1
    stratum["folded_sfs"][folded] += 1
    if derived is not None:
        stratum["oriented_sites"] += 1
        stratum["unfolded_sfs"][derived] += 1


def _serialize_stratum(stratum: dict) -> dict:
    return {
        "biallelic_sites": stratum["biallelic_sites"],
        "oriented_sites": stratum["oriented_sites"],
        "folded_sfs": {str(k): v for k, v in sorted(stratum["folded_sfs"].items())},
        "unfolded_sfs": {str(k): v for k, v in sorted(stratum["unfolded_sfs"].items())},
    }


def audit_growth_resources(vcf_path: Path, bed_path: Path, map_path: Path) -> dict:
    coding = load_coding_intervals(bed_path)
    coding_starts = {
        contig: [interval[0] for interval in intervals]
        for contig, intervals in coding.items()
    }
    recombination = load_recombination_map(map_path)
    vcf = pysam.VariantFile(str(vcf_path))
    samples = list(vcf.header.samples)
    contig_lengths = {
        contig: int(vcf.header.contigs[contig].length)
        for contig in vcf.header.contigs
        if vcf.header.contigs[contig].length is not None
    }
    if set(contig_lengths) != set(coding) or set(contig_lengths) != set(recombination):
        raise ValueError(
            "VCF, coding BED, and recombination map must contain identical contigs"
        )
    for contig, intervals in coding.items():
        if intervals and intervals[-1][1] > contig_lengths[contig]:
            raise ValueError(f"coding interval exceeds {contig} length")

    positive_rates = np.concatenate(
        [rates[rates > 0] for _, rates in recombination.values()]
    )
    rate_quartiles = (
        np.quantile(positive_rates, (0.25, 0.5, 0.75))
        if positive_rates.size
        else np.zeros(3)
    )

    def rate_label(rate: float | None) -> str:
        if rate is None:
            return "outside_map"
        if rate == 0:
            return "zero"
        return (
            f"positive_q{int(np.searchsorted(rate_quartiles, rate, side='left')) + 1}"
        )

    physical_distance_bp, physical_rate_bp, physical_joint_bp = (
        physical_stratum_lengths(contig_lengths, coding, recombination, rate_label)
    )

    distance_strata: dict[str, dict] = collections.defaultdict(_empty_stratum)
    rate_strata: dict[str, dict] = collections.defaultdict(_empty_stratum)
    joint_strata: dict[str, dict] = collections.defaultdict(_empty_stratum)
    contig_strata: dict[str, dict] = collections.defaultdict(_empty_stratum)
    windows: dict[tuple[str, int], dict] = collections.defaultdict(
        lambda: {
            "biallelic_sites": 0,
            "oriented_sites": 0,
            "pairwise_difference_sum": 0.0,
            "singleton_sites": 0,
            "recombination_rate_sum": 0.0,
            "recombination_rate_sites": 0,
        }
    )
    ancestral = collections.Counter()
    called_haplotypes = collections.Counter()
    records = biallelic = segregating = multiallelic = 0
    nonsegregating = collections.Counter()
    for record in vcf:
        validate_position(record, contig_lengths)
        records += 1
        if len(record.alleles or ()) != 2 or any(
            len(allele) != 1 for allele in record.alleles or ()
        ):
            multiallelic += int(len(record.alleles or ()) > 2)
            continue
        biallelic += 1
        alleles = diploid_alleles(record, samples, require_complete=True)
        called = len(alleles)
        alternate = int(sum(alleles))
        if called == 0:
            nonsegregating["no_called_haplotypes"] += 1
            continue
        if alternate == 0:
            nonsegregating["fixed_reference"] += 1
            continue
        if alternate == called:
            nonsegregating["fixed_alternate"] += 1
            continue
        segregating += 1
        called_haplotypes[called] += 1
        aa = record.info.get("AA")
        if isinstance(aa, tuple):
            aa = aa[0] if aa else None
        aa = str(aa).upper() if aa is not None else None
        reference, alternate_allele = (allele.upper() for allele in record.alleles)
        if aa == reference:
            derived = alternate
            ancestral["reference"] += 1
        elif aa == alternate_allele:
            derived = called - alternate
            ancestral["alternate"] += 1
        else:
            derived = None
            ancestral["missing_or_unmatched"] += 1
        folded = min(alternate, called - alternate)
        position0 = record.pos - 1
        distance = coding_distance(
            coding[record.contig], position0, coding_starts[record.contig]
        )
        rate = local_recombination_rate(recombination.get(record.contig), position0)
        dlabel = distance_label(distance)
        rlabel = rate_label(rate)
        for stratum in (
            distance_strata[dlabel],
            rate_strata[rlabel],
            joint_strata[f"{dlabel}|{rlabel}"],
            contig_strata[record.contig],
        ):
            _add_site(stratum, folded, derived)
        window = windows[(record.contig, position0 // WINDOW_SIZE)]
        window["biallelic_sites"] += 1
        window["oriented_sites"] += int(derived is not None)
        window["pairwise_difference_sum"] += (
            2 * alternate * (called - alternate) / (called * (called - 1))
        )
        window["singleton_sites"] += int(folded == 1)
        if rate is not None:
            window["recombination_rate_sum"] += rate
            window["recombination_rate_sites"] += 1

    bed_summary = {}
    map_summary = {}
    window_rows = []
    for contig in contig_lengths:
        intervals = coding[contig]
        positions, rates = recombination[contig]
        union_bp = sum(end - start for start, end in intervals)
        bed_summary[contig] = {
            "merged_intervals": len(intervals),
            "union_bp": union_bp,
            "fraction_of_contig": union_bp / contig_lengths[contig],
        }
        map_summary[contig] = {
            "points": len(positions),
            "first_position": int(positions[0]),
            "last_position": int(positions[-1]),
            "zero_rate_fraction": float(np.mean(rates == 0)),
            "rate_cm_per_mb_quantiles": dict(
                zip(
                    ("0", "0.25", "0.5", "0.75", "1"),
                    np.quantile(rates, (0, 0.25, 0.5, 0.75, 1)).tolist(),
                )
            ),
        }
        for index in range((contig_lengths[contig] + WINDOW_SIZE - 1) // WINDOW_SIZE):
            start = index * WINDOW_SIZE
            end = min((index + 1) * WINDOW_SIZE, contig_lengths[contig])
            values = windows[(contig, index)]
            span = end - start
            window_rows.append(
                {
                    "contig": contig,
                    "start": start,
                    "end": end,
                    "biallelic_sites": values["biallelic_sites"],
                    "oriented_sites": values["oriented_sites"],
                    "diversity_per_bp": values["pairwise_difference_sum"] / span,
                    "singleton_fraction": values["singleton_sites"]
                    / max(values["biallelic_sites"], 1),
                    "coding_fraction": overlap_length(intervals, start, end) / span,
                    "mean_variant_site_recombination_cm_per_mb": values[
                        "recombination_rate_sum"
                    ]
                    / max(values["recombination_rate_sites"], 1),
                }
            )

    return {
        "files": {
            "vcf": {"path": str(vcf_path), "sha256": sha256(vcf_path)},
            "coding_bed": {"path": str(bed_path), "sha256": sha256(bed_path)},
            "recombination_map": {"path": str(map_path), "sha256": sha256(map_path)},
        },
        "samples": len(samples),
        "expected_haplotypes": 2 * len(samples),
        "contig_lengths": contig_lengths,
        "sequence_length": sum(contig_lengths.values()),
        "records": records,
        "biallelic_sites": biallelic,
        "segregating_biallelic_sites": segregating,
        "nonsegregating_biallelic_sites": dict(nonsegregating),
        "multiallelic_sites": multiallelic,
        "ancestral_orientation": dict(ancestral),
        "called_haplotype_histogram": {
            str(k): v for k, v in sorted(called_haplotypes.items())
        },
        "vcf_info_fields": sorted(vcf.header.info.keys()),
        "vcf_format_fields": sorted(vcf.header.formats.keys()),
        "nonsynonymous_annotation_present": any(
            "non" in field.lower() and "syn" in field.lower()
            for field in vcf.header.info.keys()
        ),
        "positive_map_rate_quartiles_cm_per_mb": rate_quartiles.tolist(),
        "physical_bp_by_distance_stratum": physical_distance_bp,
        "physical_bp_by_recombination_stratum": physical_rate_bp,
        "physical_bp_by_joint_stratum": physical_joint_bp,
        "coding_bed": bed_summary,
        "recombination_map": map_summary,
        "distance_strata": {
            key: _serialize_stratum(value)
            for key, value in sorted(distance_strata.items())
        },
        "recombination_strata": {
            key: _serialize_stratum(value) for key, value in sorted(rate_strata.items())
        },
        "joint_strata": {
            key: _serialize_stratum(value)
            for key, value in sorted(joint_strata.items())
        },
        "contig_strata": {
            key: _serialize_stratum(value)
            for key, value in sorted(contig_strata.items())
        },
        "windows_100kb": window_rows,
    }
