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
    variant_types: collections.Counter[str] = collections.Counter()
    filters: collections.Counter[str] = collections.Counter()
    ploidies: collections.Counter[int] = collections.Counter()
    positions_by_contig: dict[str, list[int]] = collections.defaultdict(list)
    allele_counts: list[int] = []
    called_counts: list[int] = []
    phased = unphased = 0
    total = biallelic_snps = 0
    for rec in vcf:
        total += 1
        positions_by_contig[rec.contig].append(rec.pos)
        filters.update(rec.filter.keys() or ["PASS"])
        if len(rec.alleles or ()) == 2 and all(len(a) == 1 for a in rec.alleles):
            variant_types["biallelic_snp"] += 1
            biallelic_snps += 1
        elif len(rec.alleles or ()) > 2:
            variant_types["multiallelic"] += 1
        elif any(len(a) != 1 for a in rec.alleles or ()):
            variant_types["indel_or_complex"] += 1
        ac = called = 0
        for i, sample in enumerate(samples):
            call = rec.samples[sample]
            gt = call.get("GT")
            if gt is None:
                sample_missing[i] += 1
                continue
            vals = [a for a in gt if a is not None and a >= 0]
            ploidies[len(gt)] += 1
            if len(vals) != len(gt):
                sample_missing[i] += 1
            else:
                sample_called[i] += 1
            called += len(vals)
            ac += sum(a > 0 for a in vals)
            if call.phased:
                phased += 1
            else:
                unphased += 1
        if 0 < ac < called:
            allele_counts.append(ac)
            called_counts.append(called)
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
    folded_frequency = np.minimum(allele_counts, np.asarray(called_counts) - allele_counts) / np.maximum(called_counts, 1)
    return {
        "file": path.name,
        "samples": len(samples),
        "sample_names_preview": samples[:10],
        "records": total,
        "biallelic_snp_fraction": biallelic_snps / max(total, 1),
        "variant_types": dict(variant_types),
        "filters": dict(filters),
        "ploidy_call_counts": {str(k): v for k, v in ploidies.items()},
        "phased_call_fraction": phased / max(phased + unphased, 1),
        "sample_missingness_quantiles": quantiles(missing_fraction.tolist()),
        "folded_frequency_quantiles": quantiles(folded_frequency.tolist()),
        "intervariant_spacing_quantiles": quantiles(all_spacings),
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

