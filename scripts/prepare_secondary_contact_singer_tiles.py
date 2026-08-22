from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--out-root", type=Path, required=True)
    parser.add_argument("--tile-length", type=int, default=100_000)
    parser.add_argument("--sequence-length", type=int, default=100_000_000)
    args = parser.parse_args()
    if args.tile_length <= 0 or args.sequence_length % args.tile_length:
        raise ValueError("tile length must divide sequence length")
    args.out_root.mkdir(parents=True, exist_ok=True)
    samples = []
    header = ""
    current_tile = None
    destination = None
    counts: dict[int, int] = {}
    skipped: dict[int, int] = {}
    try:
        with gzip.open(args.vcf, "rt") as source:
            for line in source:
                if line.startswith("#CHROM"):
                    samples = line.rstrip().split("\t")[9:]
                    expected = [f"mainland_{index}" for index in range(1, 44)] + [
                        f"island_{index}" for index in range(1, 25)
                    ]
                    if samples != expected:
                        raise ValueError("unexpected registered sample order")
                    header = (
                        "##fileformat=VCFv4.2\n"
                        f"##contig=<ID=1,length={args.tile_length}>\n"
                        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                        + line
                    )
                    continue
                if line.startswith("#"):
                    continue
                fields = line.rstrip().split("\t")
                position = int(fields[1]) - 1
                if position >= args.sequence_length:
                    break
                tile = position // args.tile_length
                if len(fields[3]) != 1 or len(fields[4]) != 1 or "," in fields[4]:
                    skipped[tile] = skipped.get(tile, 0) + 1
                    continue
                if tile != current_tile:
                    if destination is not None:
                        destination.close()
                    tile_root = args.out_root / f"tile_{tile:04d}"
                    tile_root.mkdir(parents=True, exist_ok=True)
                    destination = (tile_root / "observed.vcf").open("w")
                    destination.write(header)
                    current_tile = tile
                fields[0] = "1"
                fields[1] = str(position % args.tile_length + 1)
                fields[2] = "."
                fields[5] = "."
                fields[6] = "PASS"
                fields[7] = "."
                fields[8] = "GT"
                destination.write("\t".join(fields) + "\n")
                counts[tile] = counts.get(tile, 0) + 1
    finally:
        if destination is not None:
            destination.close()
    expected_tiles = args.sequence_length // args.tile_length
    missing = sorted(set(range(expected_tiles)) - set(counts))
    if len(samples) != 67 or missing:
        raise ValueError(f"invalid tiling: samples={len(samples)}, empty_tiles={missing[:10]}")
    report = {
        "source": str(args.vcf),
        "output_root": str(args.out_root),
        "tile_length": args.tile_length,
        "sequence_length": args.sequence_length,
        "tiles": expected_tiles,
        "total_records": sum(counts.values()),
        "total_skipped_non_biallelic": sum(skipped.values()),
        "minimum_records_per_tile": min(counts.values()),
        "maximum_records_per_tile": max(counts.values()),
        "records_per_tile": {str(key): value for key, value in sorted(counts.items())},
        "skipped_per_tile": {str(key): value for key, value in sorted(skipped.items())},
    }
    (args.out_root / "manifest.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({key: value for key, value in report.items() if not key.endswith("per_tile")}, indent=2))


if __name__ == "__main__":
    main()
