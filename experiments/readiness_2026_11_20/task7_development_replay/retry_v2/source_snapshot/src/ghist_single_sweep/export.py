"""Canonical biallelic VCF export for genotype-only sweep statistics.

SLiM mutation IDs are not nucleotide alleles. Relabel binary states A/C for both
our predictor and RAiSD, preserving identical genotypes, sites and exposure.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np


def write_biallelic_vcf(ts, path: Path, chromosome: str = "7") -> dict:
    length = int(ts.sequence_length)
    if length != ts.sequence_length or length <= 0:
        raise ValueError("Sweep VCF export requires an integer sequence length")
    samples = list(ts.samples())
    if len(samples) % 2:
        raise ValueError("Sweep VCF export requires diploid samples")
    individuals = [ts.node(node).individual for node in samples]
    if any(i >= 0 for i in individuals):
        groups = {}
        for node, individual in zip(samples, individuals, strict=True):
            groups.setdefault(individual, []).append(node)
        if -1 in groups or any(len(nodes) != 2 for nodes in groups.values()):
            raise ValueError("Each sampled individual must contain two haplotypes")
        samples = [node for nodes in groups.values() for node in nodes]
    accepted = skipped = non_nucleotide = 0
    previous = -1
    with path.open("w") as handle:
        handle.write(
            f"##fileformat=VCFv4.2\n##contig=<ID={chromosome},length={length}>\n"
        )
        handle.write(
            '##source=GHIST_binary_genotype_export\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        )
        handle.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
            + "\t".join(f"sample_{i}" for i in range(len(samples) // 2))
            + "\n"
        )
        for variant in ts.variants(samples=samples):
            calls = variant.genotypes
            if (
                len(variant.alleles) != 2
                or not np.all((calls == 0) | (calls == 1))
                or int(calls.sum()) in (0, len(calls))
            ):
                skipped += 1
                continue
            position = int(variant.site.position)
            if position != variant.site.position or not previous < position < length:
                raise ValueError(
                    "Variants must have unique integer coordinates within the contig"
                )
            previous = position
            non_nucleotide += any(
                a not in {"A", "C", "G", "T"} for a in variant.alleles
            )
            pairs = (f"{calls[i]}|{calls[i + 1]}" for i in range(0, len(calls), 2))
            handle.write(
                f"{chromosome}\t{position + 1}\t.\tA\tC\t.\tPASS\t.\tGT\t"
                + "\t".join(pairs)
                + "\n"
            )
            accepted += 1
    return {
        "contig_length": length,
        "diploid_samples": len(samples) // 2,
        "accepted_biallelic_polymorphic_sites": accepted,
        "skipped_non_biallelic_or_fixed_sites": skipped,
        "accepted_sites_with_non_nucleotide_source_labels": non_nucleotide,
        "allele_encoding": "binary states relabeled A/C identically for predictor and RAiSD",
        "coordinate_system": "VCF site = tree site + 1; contig length unchanged",
    }
