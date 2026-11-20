from __future__ import annotations

from pathlib import Path
from typing import Collection

import numpy as np
import pysam

from .audit import coding_distance, distance_label, load_coding_intervals


def distance_stratum_segments(
    intervals: list[tuple[int, int]],
    length: int,
    included_labels: Collection[str],
) -> list[tuple[int, int]]:
    """Return half-open genome segments belonging to selected distance bands."""
    breakpoints = {0, length}
    for start, end in intervals:
        breakpoints.update((start, end))
        for threshold in (1_000, 10_000, 100_000):
            breakpoints.add(max(0, start - threshold))
            breakpoints.add(min(length, end + threshold))
    starts = [interval[0] for interval in intervals]
    selected = []
    ordered = sorted(breakpoints)
    for left, right in zip(ordered[:-1], ordered[1:]):
        if right <= left:
            continue
        if distance_label(coding_distance(intervals, left, starts)) in included_labels:
            if selected and left == selected[-1][1]:
                selected[-1] = (selected[-1][0], right)
            else:
                selected.append((left, right))
    return selected


def _overlap(intervals: list[tuple[int, int]], start: int, end: int) -> int:
    total = 0
    for left, right in intervals:
        if left >= end:
            break
        if right > start:
            total += max(0, min(right, end) - max(left, start))
    return total


def extract_block_spectra(
    vcf_path: Path,
    bed_path: Path,
    *,
    block_size: int = 5_000_000,
    included_labels: Collection[str] = (
        "noncoding_1kb_10kb",
        "noncoding_10kb_100kb",
    ),
) -> dict[str, np.ndarray]:
    """Extract unfolded SFS and exact exposure for fixed physical blocks."""
    if block_size <= 0:
        raise ValueError("block size must be positive")
    coding = load_coding_intervals(bed_path)
    coding_starts = {
        contig: [interval[0] for interval in intervals]
        for contig, intervals in coding.items()
    }
    vcf = pysam.VariantFile(str(vcf_path))
    samples = list(vcf.header.samples)
    haplotypes = 2 * len(samples)
    contig_lengths = {
        contig: int(vcf.header.contigs[contig].length)
        for contig in vcf.header.contigs
        if vcf.header.contigs[contig].length is not None
    }
    rows = []
    block_lookup = {}
    for contig, length in contig_lengths.items():
        selected = distance_stratum_segments(
            coding[contig], length, included_labels
        )
        for start in range(0, length, block_size):
            end = min(start + block_size, length)
            index = len(rows)
            block_lookup[(contig, start // block_size)] = index
            rows.append(
                {
                    "contig": contig,
                    "start": start,
                    "end": end,
                    "exposure": _overlap(selected, start, end),
                    "spectrum": np.zeros(haplotypes + 1, dtype=np.int64),
                }
            )

    accepted = 0
    for record in vcf:
        if len(record.alleles or ()) != 2 or any(
            len(allele) != 1 for allele in record.alleles or ()
        ):
            continue
        alleles = []
        for sample in samples:
            genotype = record.samples[sample].get("GT")
            if genotype is not None:
                alleles.extend(
                    allele for allele in genotype if allele is not None and allele in (0, 1)
                )
        if len(alleles) != haplotypes:
            continue
        alternate_count = int(sum(alleles))
        if alternate_count in (0, haplotypes):
            continue
        position0 = record.pos - 1
        label = distance_label(
            coding_distance(
                coding[record.contig], position0, coding_starts[record.contig]
            )
        )
        if label not in included_labels:
            continue
        ancestral = record.info.get("AA")
        if isinstance(ancestral, tuple):
            ancestral = ancestral[0] if ancestral else None
        ancestral = str(ancestral).upper() if ancestral is not None else None
        reference, alternate = (allele.upper() for allele in record.alleles)
        if ancestral == reference:
            derived = alternate_count
        elif ancestral == alternate:
            derived = haplotypes - alternate_count
        else:
            continue
        block = block_lookup[(record.contig, position0 // block_size)]
        rows[block]["spectrum"][derived] += 1
        accepted += 1
    vcf.close()
    return {
        "spectra": np.stack([row["spectrum"] for row in rows]),
        "exposure_bp": np.asarray([row["exposure"] for row in rows], dtype=np.int64),
        "contig": np.asarray([row["contig"] for row in rows]),
        "start": np.asarray([row["start"] for row in rows], dtype=np.int64),
        "end": np.asarray([row["end"] for row in rows], dtype=np.int64),
        "block_size": np.asarray(block_size, dtype=np.int64),
        "haplotypes": np.asarray(haplotypes, dtype=np.int64),
        "accepted_sites": np.asarray(accepted, dtype=np.int64),
        "included_labels": np.asarray(sorted(included_labels)),
    }
