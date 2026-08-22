from __future__ import annotations

import argparse
import collections
import hashlib
import json
from pathlib import Path

import numpy as np
import pysam

from audit_vcf import audit


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def population_audit(path: Path) -> dict:
    vcf = pysam.VariantFile(str(path))
    samples = list(vcf.header.samples)
    groups = {
        "mainland": [sample for sample in samples if sample.startswith("mainland_")],
        "island": [sample for sample in samples if sample.startswith("island_")],
    }
    if set(groups["mainland"]) | set(groups["island"]) != set(samples):
        raise ValueError("sample names do not partition into mainland and island")
    statistics = collections.Counter()
    joint_counts: collections.Counter[tuple[int, int]] = collections.Counter()
    pi_sum = np.zeros(2, dtype=float)
    dxy_sum = fst_numerator = 0.0
    windows = collections.defaultdict(lambda: np.zeros(6, dtype=float))
    called_histograms = [collections.Counter(), collections.Counter()]
    for record in vcf:
        if len(record.alleles or ()) != 2 or any(len(a) != 1 for a in record.alleles):
            continue
        alt_counts = []
        called_counts = []
        all_phased = True
        for group in groups.values():
            alleles = []
            for sample in group:
                call = record.samples[sample]
                genotype = call.get("GT")
                if genotype is not None:
                    alleles.extend(a for a in genotype if a is not None and a >= 0)
                    all_phased &= bool(call.phased)
            alt_counts.append(sum(a == 1 for a in alleles))
            called_counts.append(len(alleles))
        ancestral = record.info.get("AA")
        if isinstance(ancestral, tuple):
            ancestral = ancestral[0]
        if ancestral == record.ref:
            statistics["ancestral_reference"] += 1
            derived = alt_counts
        elif ancestral == record.alts[0]:
            statistics["ancestral_alternate"] += 1
            derived = [n - count for n, count in zip(called_counts, alt_counts)]
        else:
            statistics["ancestral_invalid"] += 1
            continue
        statistics["biallelic_snps"] += 1
        statistics["fully_phased_snps"] += int(all_phased)
        for index, called in enumerate(called_counts):
            called_histograms[index][called] += 1
        if 0 < derived[0] < called_counts[0]:
            statistics["segregating_mainland"] += 1
        if 0 < derived[1] < called_counts[1]:
            statistics["segregating_island"] += 1
        if 0 < derived[0] < called_counts[0] and derived[1] in (0, called_counts[1]):
            statistics["private_mainland"] += 1
        if 0 < derived[1] < called_counts[1] and derived[0] in (0, called_counts[0]):
            statistics["private_island"] += 1
        if all(0 < count < called for count, called in zip(derived, called_counts)):
            statistics["shared_polymorphism"] += 1
        if (derived[0] == 0 and derived[1] == called_counts[1]) or (
            derived[1] == 0 and derived[0] == called_counts[0]
        ):
            statistics["fixed_difference"] += 1
        if all(called > 1 for called in called_counts):
            frequencies = np.asarray(derived) / np.asarray(called_counts)
            pi = 2 * frequencies * (1 - frequencies) * np.asarray(called_counts) / (
                np.asarray(called_counts) - 1
            )
            dxy = frequencies[0] * (1 - frequencies[1]) + frequencies[1] * (
                1 - frequencies[0]
            )
            pi_sum += pi
            dxy_sum += dxy
            fst_numerator += dxy - 0.5 * pi.sum()
            window = windows[(record.pos - 1) // 1_000_000]
            window[:2] += pi
            window[2] += dxy
            window[3] += int(0 < derived[0] < called_counts[0])
            window[4] += int(0 < derived[1] < called_counts[1])
            window[5] += 1
        if called_counts == [86, 48]:
            joint_counts[(derived[0], derived[1])] += 1
    denominator = max(statistics["biallelic_snps"], 1)
    return {
        "file_sha256": sha256(path),
        "populations": {name: members for name, members in groups.items()},
        "diploid_sample_counts": {name: len(members) for name, members in groups.items()},
        "haploid_sample_counts": {name: 2 * len(members) for name, members in groups.items()},
        "statistics": dict(statistics),
        "fully_phased_fraction": statistics["fully_phased_snps"] / denominator,
        "mean_pairwise_diversity_per_base": (pi_sum / 100_000_000).tolist(),
        "mean_dxy_per_base": dxy_sum / 100_000_000,
        "hudson_fst": fst_numerator / max(dxy_sum, 1e-12),
        "called_haplotype_histograms": [
            {str(key): value for key, value in sorted(histogram.items())}
            for histogram in called_histograms
        ],
        "complete_joint_sfs_nonzero_cells": len(joint_counts),
        "one_mb_windows": {
            str(key): value.tolist() for key, value in sorted(windows.items())
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = {"generic": audit(args.vcf), "population": population_audit(args.vcf)}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
