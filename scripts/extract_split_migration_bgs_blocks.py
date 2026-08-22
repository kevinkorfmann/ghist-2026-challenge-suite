from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pysam

from ghist_growth_bgs.audit import coding_distance, load_coding_intervals
from ghist_growth_bgs.blocks import distance_stratum_segments


MIDDLE_LABELS = {"noncoding_1kb_10kb", "noncoding_10kb_100kb"}


def overlap_length(
    intervals: list[tuple[int, int]], start: int, end: int
) -> int:
    return sum(
        max(0, min(right, end) - max(left, start))
        for left, right in intervals
        if left < end and right > start
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--coding-bed", type=Path, required=True)
    parser.add_argument("--block-size", type=int, default=1_000_000)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args()
    coding = load_coding_intervals(args.coding_bed)
    starts = {contig: [left for left, _ in values] for contig, values in coding.items()}
    vcf = pysam.VariantFile(str(args.vcf))
    afr = [sample for sample in vcf.header.samples if sample.startswith("AFR_")]
    eur = [sample for sample in vcf.header.samples if sample.startswith("EUR_")]
    shape = (2 * len(afr) + 1, 2 * len(eur) + 1)
    contig_lengths = {
        name: int(vcf.header.contigs[name].length) for name in vcf.header.contigs
    }
    rows = []
    lookup = {}
    for contig, length in contig_lengths.items():
        selected = distance_stratum_segments(coding[contig], length, MIDDLE_LABELS)
        for index in range(length // args.block_size):
            left = index * args.block_size
            right = left + args.block_size
            lookup[(contig, index)] = len(rows)
            rows.append(
                {
                    "contig": contig,
                    "block_index": index,
                    "left": left,
                    "right": right,
                    "exposure_bp": overlap_length(selected, left, right),
                }
            )
    spectra = np.zeros((len(rows), *shape), dtype=np.int32)
    accepted = skipped = 0
    for record in vcf:
        if len(record.alleles or ()) != 2 or any(len(a) != 1 for a in record.alleles):
            skipped += 1
            continue
        position = record.pos - 1
        row_index = lookup.get((record.contig, position // args.block_size))
        if row_index is None:
            continue
        distance = coding_distance(coding[record.contig], position, starts[record.contig])
        if not (distance is not None and 1_000 < distance <= 100_000):
            continue
        counts = []
        called = []
        for samples in (afr, eur):
            alleles = [
                allele
                for sample in samples
                for allele in (record.samples[sample].get("GT") or ())
                if allele in (0, 1)
            ]
            counts.append(sum(alleles))
            called.append(len(alleles))
        if tuple(called) != (shape[0] - 1, shape[1] - 1):
            skipped += 1
            continue
        ancestral = record.info.get("AA")
        if isinstance(ancestral, tuple):
            ancestral = ancestral[0]
        if ancestral == record.alts[0]:
            counts = [n - k for n, k in zip(called, counts)]
        elif ancestral != record.ref:
            skipped += 1
            continue
        spectra[row_index, counts[0], counts[1]] += 1
        accepted += 1
    args.out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.out,
        spectra=spectra,
        exposure_bp=np.asarray([row["exposure_bp"] for row in rows], dtype=np.int64),
        contig=np.asarray([row["contig"] for row in rows]),
        block_index=np.asarray([row["block_index"] for row in rows], dtype=np.int64),
    )
    payload = {
        "state": "COMPLETED",
        "block_size": args.block_size,
        "blocks": len(rows),
        "blocks_by_contig": {
            contig: sum(row["contig"] == contig for row in rows)
            for contig in contig_lengths
        },
        "physical_exposure_bp": sum(row["exposure_bp"] for row in rows),
        "accepted_sites": accepted,
        "skipped_records": skipped,
        "terminal_partial_blocks_excluded": True,
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()

