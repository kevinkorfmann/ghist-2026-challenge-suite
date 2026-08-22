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


def scan_secondary_contact_vcf(path: str | Path) -> tuple[np.ndarray, dict]:
    """Scan the registered GT-only VCF with strict layout validation.

    The challenge stores 43 mainland samples followed by 24 island samples,
    with complete diploid genotypes. Avoiding millions of per-sample call
    objects makes this exact scanner substantially faster than generic parsing.
    """
    path = Path(path)
    spectrum = np.zeros(
        (POPULATION1_HAPLOTYPES + 1, POPULATION2_HAPLOTYPES + 1), dtype=np.int64
    )
    totals = np.zeros(12, dtype=np.float64)
    windows = np.zeros((100, 7), dtype=np.float64)
    records = ancestral_reference = ancestral_alternate = 0
    invalid_ancestral = layout_errors = skipped_non_biallelic = 0
    nonmonotone = duplicate = 0
    first_position = last_position = previous_position = None
    valid_genotype_bytes = set(b"01|\t")
    with gzip.open(path, "rb") as handle:
        for line in handle:
            if line.startswith(b"#"):
                continue
            fields = line.rstrip(b"\r\n").split(b"\t", 9)
            if len(fields) != 10 or fields[8] != b"GT":
                layout_errors += 1
                continue
            genotypes = fields[9]
            if (
                len(genotypes) != 4 * (POPULATION1_DIPLOIDS + POPULATION2_DIPLOIDS) - 1
                or any(value not in valid_genotype_bytes for value in genotypes)
                or genotypes[POPULATION1_BYTES] != 9
            ):
                layout_errors += 1
                continue
            position = int(fields[1])
            if first_position is None:
                first_position = position
            if previous_position is not None:
                nonmonotone += int(position < previous_position)
                duplicate += int(position == previous_position)
            previous_position = last_position = position
            ref, alt, info = fields[3], fields[4], fields[7]
            if len(ref) != 1 or len(alt) != 1 or b"," in alt:
                skipped_non_biallelic += 1
                continue
            marker = info.find(b"AA=")
            ancestral = (
                info[marker + 3 :].split(b";", 1)[0] if marker >= 0 else b""
            )
            count1 = genotypes[:POPULATION1_BYTES].count(b"1")
            count2 = genotypes[POPULATION1_BYTES + 1 :].count(b"1")
            if ancestral == ref:
                ancestral_reference += 1
            elif ancestral == alt:
                ancestral_alternate += 1
                count1 = POPULATION1_HAPLOTYPES - count1
                count2 = POPULATION2_HAPLOTYPES - count2
            else:
                invalid_ancestral += 1
                continue
            spectrum[count1, count2] += 1
            records += 1
            p1 = count1 / POPULATION1_HAPLOTYPES
            p2 = count2 / POPULATION2_HAPLOTYPES
            pi1 = 2 * p1 * (1 - p1) * POPULATION1_HAPLOTYPES / (
                POPULATION1_HAPLOTYPES - 1
            )
            pi2 = 2 * p2 * (1 - p2) * POPULATION2_HAPLOTYPES / (
                POPULATION2_HAPLOTYPES - 1
            )
            dxy = p1 * (1 - p2) + p2 * (1 - p1)
            flags = (
                int(0 < count1 < POPULATION1_HAPLOTYPES),
                int(0 < count2 < POPULATION2_HAPLOTYPES),
                int(0 < count1 < POPULATION1_HAPLOTYPES and count2 in (0, POPULATION2_HAPLOTYPES)),
                int(0 < count2 < POPULATION2_HAPLOTYPES and count1 in (0, POPULATION1_HAPLOTYPES)),
                int(0 < count1 < POPULATION1_HAPLOTYPES and 0 < count2 < POPULATION2_HAPLOTYPES),
                int(
                    (count1 == 0 and count2 == POPULATION2_HAPLOTYPES)
                    or (count2 == 0 and count1 == POPULATION1_HAPLOTYPES)
                ),
            )
            totals[:3] += (pi1, pi2, dxy)
            totals[3] += dxy - 0.5 * (pi1 + pi2)
            totals[4:10] += flags
            totals[10] += abs(p1 - p2)
            totals[11] += 1
            window = min((position - 1) // 1_000_000, 99)
            windows[window] += (pi1, pi2, dxy, *flags[:3], 1)
    if records == 0 or layout_errors:
        raise ValueError(
            f"VCF layout validation failed: records={records}, layout_errors={layout_errors}"
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
        "diploid_sample_counts": {"mainland": 43, "island": 24},
        "haploid_sample_counts": {"mainland": 86, "island": 48},
        "fully_called_fraction": 1.0,
        "fully_phased_fraction": 1.0,
        "mean_pairwise_diversity_per_base": [
            totals[0] / 100_000_000,
            totals[1] / 100_000_000,
        ],
        "mean_dxy_per_base": totals[2] / 100_000_000,
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
