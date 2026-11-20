from __future__ import annotations

import gzip
import hashlib
from pathlib import Path

import numpy as np


POPULATION1_DIPLOIDS = 43
POPULATION2_DIPLOIDS = 24
POPULATION1_HAPLOTYPES = 2 * POPULATION1_DIPLOIDS
POPULATION2_HAPLOTYPES = 2 * POPULATION2_DIPLOIDS
POPULATION1_BYTES = 4 * POPULATION1_DIPLOIDS - 1


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def scan_secondary_contact_vcf(
    path: str | Path, *, sequence_length: int = 100_000_000
) -> tuple[np.ndarray, dict]:
    """Scan complete GT-only VCFs using population membership from the header.

    Sample names must start with mainland_ or island_. Both population sizes
    and column order are inferred; complete phased or unphased diploids work.
    """
    path = Path(path)
    if not isinstance(sequence_length, int) or sequence_length <= 0:
        raise ValueError("Sequence length must be a positive integer")
    spectrum = None
    sample_names = None
    n1 = n2 = h1 = h2 = total_samples = split_bytes = 0
    contiguous = False
    indices1, indices2 = [], []
    totals = np.zeros(12, dtype=np.float64)
    windows = np.zeros(((sequence_length + 999_999) // 1_000_000, 7), dtype=np.float64)
    records = ancestral_reference = ancestral_alternate = 0
    invalid_ancestral = layout_errors = skipped_non_biallelic = 0
    nonmonotone = duplicate = 0
    first_position = last_position = previous_position = None
    valid_genotype_bytes = set(b"01|/\t")
    phased_genotypes = 0
    contig = None
    opener = gzip.open if path.suffix in {".gz", ".bgz"} else open
    with opener(path, "rb") as handle:
        for line in handle:
            if line.startswith(b"#CHROM\t"):
                if sample_names is not None:
                    raise ValueError("Duplicate VCF sample header")
                sample_names = line.rstrip().split(b"\t")[9:]
                indices1 = [
                    i
                    for i, name in enumerate(sample_names)
                    if name.startswith(b"mainland_")
                ]
                indices2 = [
                    i
                    for i, name in enumerate(sample_names)
                    if name.startswith(b"island_")
                ]
                n1, n2 = len(indices1), len(indices2)
                total_samples = n1 + n2
                if (
                    not n1
                    or not n2
                    or total_samples != len(sample_names)
                    or len(set(sample_names)) != total_samples
                ):
                    raise ValueError(
                        "Header requires unique mainland_ and island_ sample names"
                    )
                h1, h2, split_bytes = 2 * n1, 2 * n2, 4 * n1 - 1
                spectrum = np.zeros((h1 + 1, h2 + 1), dtype=np.int64)
                contiguous = indices1 == list(range(n1))
                continue
            if line.startswith(b"#"):
                continue
            if sample_names is None:
                raise ValueError("VCF sample header is missing")
            fields = line.rstrip(b"\r\n").split(b"\t", 9)
            if len(fields) != 10 or fields[8] != b"GT":
                layout_errors += 1
                continue
            ref, alt, info = fields[3], fields[4], fields[7]
            if len(ref) != 1 or len(alt) != 1 or b"," in alt:
                skipped_non_biallelic += 1
                continue
            genotypes = fields[9]
            if (
                len(genotypes) != 4 * (n1 + n2) - 1
                or any(value not in valid_genotype_bytes for value in genotypes)
                or any(value not in b"|/" for value in genotypes[1::4])
                or genotypes[3::4] != b"\t" * (n1 + n2 - 1)
                or any(value not in b"01" for value in genotypes[0::4])
                or any(value not in b"01" for value in genotypes[2::4])
            ):
                layout_errors += 1
                continue
            if contig is None:
                contig = fields[0]
            if fields[0] != contig:
                raise ValueError(
                    "Secondary-contact fast scanner requires a single contig"
                )
            phased_genotypes += genotypes[1::4].count(b"|")
            position = int(fields[1])
            if not 1 <= position <= sequence_length:
                raise ValueError(
                    "VCF positions must be 1-based within the specified region"
                )
            if first_position is None:
                first_position = position
            if previous_position is not None:
                nonmonotone += int(position < previous_position)
                duplicate += int(position == previous_position)
            previous_position = last_position = position
            annotations = [
                item[3:] for item in info.split(b";") if item.startswith(b"AA=")
            ]
            ancestral = annotations[0] if len(annotations) == 1 else b""
            if contiguous:
                count1 = genotypes[:split_bytes].count(b"1")
                count2 = genotypes[split_bytes + 1 :].count(b"1")
            else:
                count1 = sum(
                    (genotypes[4 * i] == 49) + (genotypes[4 * i + 2] == 49)
                    for i in indices1
                )
                count2 = sum(
                    (genotypes[4 * i] == 49) + (genotypes[4 * i + 2] == 49)
                    for i in indices2
                )
            if ancestral == ref:
                ancestral_reference += 1
            elif ancestral == alt:
                ancestral_alternate += 1
                count1 = h1 - count1
                count2 = h2 - count2
            else:
                invalid_ancestral += 1
                continue
            spectrum[count1, count2] += 1
            records += 1
            p1 = count1 / h1
            p2 = count2 / h2
            pi1 = 2 * p1 * (1 - p1) * h1 / (h1 - 1)
            pi2 = 2 * p2 * (1 - p2) * h2 / (h2 - 1)
            dxy = p1 * (1 - p2) + p2 * (1 - p1)
            flags = (
                int(0 < count1 < h1),
                int(0 < count2 < h2),
                int(0 < count1 < h1 and count2 in (0, h2)),
                int(0 < count2 < h2 and count1 in (0, h1)),
                int(0 < count1 < h1 and 0 < count2 < h2),
                int((count1 == 0 and count2 == h2) or (count2 == 0 and count1 == h1)),
            )
            totals[:3] += (pi1, pi2, dxy)
            totals[3] += dxy - 0.5 * (pi1 + pi2)
            totals[4:10] += flags
            totals[10] += abs(p1 - p2)
            totals[11] += 1
            window = (position - 1) // 1_000_000
            windows[window] += (pi1, pi2, dxy, *flags[:3], 1)
    if records == 0 or layout_errors or invalid_ancestral or nonmonotone or duplicate:
        raise ValueError(
            f"VCF layout validation failed: records={records}, layout_errors={layout_errors}, invalid_ancestral={invalid_ancestral}, nonmonotone={nonmonotone}, duplicates={duplicate}"
        )
    metadata = {
        "file_sha256": _sha256(path),
        "records": records,
        "first_position": first_position,
        "last_position": last_position,
        "duplicate_adjacent_positions": duplicate,
        "nonmonotone_positions": nonmonotone,
        "ancestral_reference_records": ancestral_reference,
        "ancestral_alternate_records": ancestral_alternate,
        "ancestral_invalid_records": invalid_ancestral,
        "skipped_non_biallelic_records": skipped_non_biallelic,
        "diploid_sample_counts": {"mainland": n1, "island": n2},
        "haploid_sample_counts": {"mainland": h1, "island": h2},
        "fully_called_fraction": 1.0,
        "fully_phased_fraction": phased_genotypes / (total_samples * records),
        "mean_pairwise_diversity_per_base": [
            totals[0] / sequence_length,
            totals[1] / sequence_length,
        ],
        "mean_dxy_per_base": totals[2] / sequence_length,
        "hudson_fst": totals[3] / max(totals[2], 1e-12),
        "segregating_mainland": int(totals[4]),
        "segregating_island": int(totals[5]),
        "private_mainland": int(totals[6]),
        "private_island": int(totals[7]),
        "shared_polymorphism": int(totals[8]),
        "fixed_difference": int(totals[9]),
        "mean_absolute_frequency_difference": totals[10] / totals[11],
        "complete_joint_sfs_nonzero_cells": int(np.count_nonzero(spectrum)),
        "one_mb_windows": windows.tolist(),
    }
    return spectrum, metadata
