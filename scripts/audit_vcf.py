from __future__ import annotations

import argparse
import collections
import json
from pathlib import Path

import numpy as np
import pysam


def quantiles(values: list[float]) -> dict[str, float | None]:
    if not values:
        return {str(q): None for q in (0, 0.01, 0.1, 0.5, 0.9, 0.99, 1)}
    array = np.asarray(values, dtype=float)
    return {str(q): float(np.quantile(array, q)) for q in (0, 0.01, 0.1, 0.5, 0.9, 0.99, 1)}


def audit(path: Path) -> dict:
    vcf = pysam.VariantFile(str(path))
    samples = list(vcf.header.samples)
    sample_called = np.zeros(len(samples), dtype=np.int64)
    sample_missing = np.zeros(len(samples), dtype=np.int64)
    sample_heterozygous = np.zeros(len(samples), dtype=np.int64)
    sample_complete_diploid = np.zeros(len(samples), dtype=np.int64)
    variant_types: collections.Counter[str] = collections.Counter()
    filters: collections.Counter[str] = collections.Counter()
    ploidies: collections.Counter[int] = collections.Counter()
    positions_by_contig: dict[str, list[int]] = collections.defaultdict(list)
    allele_counts: list[int] = []
    called_counts: list[int] = []
    alt_allele_count_histogram: collections.Counter[int] = collections.Counter()
    minor_allele_count_histogram: collections.Counter[int] = collections.Counter()
    called_haplotype_histogram: collections.Counter[int] = collections.Counter()
    genotype_classes: collections.Counter[str] = collections.Counter()
    substitution_classes: collections.Counter[str] = collections.Counter()
    window_records: collections.Counter[tuple[str, int]] = collections.Counter()
    window_snps: collections.Counter[tuple[str, int]] = collections.Counter()
    window_missing_calls: collections.Counter[tuple[str, int]] = collections.Counter()
    window_total_calls: collections.Counter[tuple[str, int]] = collections.Counter()
    quality_values: list[float] = []
    depth_values: list[float] = []
    genotype_quality_values: list[float] = []
    site_call_rates: list[float] = []
    site_heterozygosity: list[float] = []
    duplicate_positions = nonmonotone_positions = 0
    last_position: dict[str, int] = {}
    phased = unphased = 0
    total = biallelic_snps = 0
    for rec in vcf:
        total += 1
        positions_by_contig[rec.contig].append(rec.pos)
        previous = last_position.get(rec.contig)
        if previous is not None:
            duplicate_positions += int(rec.pos == previous)
            nonmonotone_positions += int(rec.pos < previous)
        last_position[rec.contig] = rec.pos
        window_key = (rec.contig, (rec.pos - 1) // 1_000_000)
        window_records[window_key] += 1
        filters.update(rec.filter.keys() or ["PASS"])
        if rec.qual is not None and np.isfinite(rec.qual):
            quality_values.append(float(rec.qual))
        if len(rec.alleles or ()) == 2 and all(len(a) == 1 for a in rec.alleles):
            variant_types["biallelic_snp"] += 1
            biallelic_snps += 1
            window_snps[window_key] += 1
            ref, alt = (a.upper() for a in rec.alleles)
            pair = "".join(sorted((ref, alt)))
            substitution_classes[
                "transition" if pair in {"AG", "CT"} else "transversion"
            ] += 1
        elif len(rec.alleles or ()) > 2:
            variant_types["multiallelic"] += 1
        elif any(len(a) != 1 for a in rec.alleles or ()):
            variant_types["indel_or_complex"] += 1
        ac = called = complete_calls = heterozygous_calls = 0
        for i, sample in enumerate(samples):
            call = rec.samples[sample]
            gt = call.get("GT")
            window_total_calls[window_key] += 1
            if gt is None:
                sample_missing[i] += 1
                genotype_classes["missing"] += 1
                window_missing_calls[window_key] += 1
                continue
            vals = [a for a in gt if a is not None and a >= 0]
            ploidies[len(gt)] += 1
            if len(vals) != len(gt):
                sample_missing[i] += 1
                genotype_classes["partial_or_missing"] += 1
                window_missing_calls[window_key] += 1
            else:
                sample_called[i] += 1
                complete_calls += 1
            called += len(vals)
            ac += sum(a > 0 for a in vals)
            if len(vals) == 2:
                sample_complete_diploid[i] += 1
                if vals[0] != vals[1]:
                    sample_heterozygous[i] += 1
                    heterozygous_calls += 1
                    genotype_classes["heterozygous"] += 1
                elif vals[0] == 0:
                    genotype_classes["homozygous_reference"] += 1
                else:
                    genotype_classes["homozygous_alternate"] += 1
            elif vals:
                genotype_classes["other_complete_ploidy"] += 1
            depth = call.get("DP")
            if depth is not None and np.isfinite(depth):
                depth_values.append(float(depth))
            genotype_quality = call.get("GQ")
            if genotype_quality is not None and np.isfinite(genotype_quality):
                genotype_quality_values.append(float(genotype_quality))
            if call.phased:
                phased += 1
            else:
                unphased += 1
        site_call_rates.append(complete_calls / max(len(samples), 1))
        site_heterozygosity.append(heterozygous_calls / max(complete_calls, 1))
        if 0 < ac < called:
            allele_counts.append(ac)
            called_counts.append(called)
            alt_allele_count_histogram[ac] += 1
            minor_allele_count_histogram[min(ac, called - ac)] += 1
            called_haplotype_histogram[called] += 1
    contigs = {}
    all_spacings = []
    for contig, positions in positions_by_contig.items():
        spacings = np.diff(np.asarray(positions, dtype=np.int64)).tolist()
        all_spacings.extend(spacings)
        header_length = vcf.header.contigs[contig].length if contig in vcf.header.contigs else None
        contigs[contig] = {
            "header_length": header_length,
            "first_position": min(positions),
            "last_position": max(positions),
            "variant_count": len(positions),
            "spacing_quantiles": quantiles(spacings),
        }
    missing_fraction = sample_missing / np.maximum(sample_called + sample_missing, 1)
    sample_heterozygosity = sample_heterozygous / np.maximum(sample_complete_diploid, 1)
    folded_frequency = np.minimum(allele_counts, np.asarray(called_counts) - allele_counts) / np.maximum(called_counts, 1)
    one_mb_windows = []
    for contig, window_index in sorted(window_records):
        key = (contig, window_index)
        one_mb_windows.append(
            {
                "contig": contig,
                "start": window_index * 1_000_000,
                "end": (window_index + 1) * 1_000_000,
                "records": window_records[key],
                "biallelic_snps": window_snps[key],
                "missing_call_fraction": window_missing_calls[key]
                / max(window_total_calls[key], 1),
            }
        )
    return {
        "file": path.name,
        "samples": len(samples),
        "sample_names_preview": samples[:10],
        "records": total,
        "biallelic_snp_fraction": biallelic_snps / max(total, 1),
        "variant_types": dict(variant_types),
        "filters": dict(filters),
        "ploidy_call_counts": {str(k): v for k, v in ploidies.items()},
        "genotype_class_counts": dict(genotype_classes),
        "phased_call_fraction": phased / max(phased + unphased, 1),
        "duplicate_adjacent_positions": duplicate_positions,
        "nonmonotone_positions": nonmonotone_positions,
        "substitution_class_counts": dict(substitution_classes),
        "site_quality_quantiles": quantiles(quality_values),
        "call_depth_quantiles": quantiles(depth_values),
        "genotype_quality_quantiles": quantiles(genotype_quality_values),
        "site_call_rate_quantiles": quantiles(site_call_rates),
        "site_heterozygosity_quantiles": quantiles(site_heterozygosity),
        "sample_missingness_quantiles": quantiles(missing_fraction.tolist()),
        "sample_heterozygosity_quantiles": quantiles(sample_heterozygosity.tolist()),
        "folded_frequency_quantiles": quantiles(folded_frequency.tolist()),
        "alt_allele_count_histogram": {
            str(k): v for k, v in sorted(alt_allele_count_histogram.items())
        },
        "minor_allele_count_histogram": {
            str(k): v for k, v in sorted(minor_allele_count_histogram.items())
        },
        "called_haplotype_histogram": {
            str(k): v for k, v in sorted(called_haplotype_histogram.items())
        },
        "intervariant_spacing_quantiles": quantiles(all_spacings),
        "one_mb_windows": one_mb_windows,
        "contigs": contigs,
        "format_fields": sorted(vcf.header.formats.keys()),
        "info_fields": sorted(vcf.header.info.keys()),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    report = audit(args.vcf)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
