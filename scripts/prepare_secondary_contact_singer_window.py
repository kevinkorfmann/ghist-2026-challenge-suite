from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--start", type=int, required=True)
    parser.add_argument("--end", type=int, required=True)
    args = parser.parse_args()
    if not 0 <= args.start < args.end:
        raise ValueError("expected 0 <= start < end")
    samples = []
    records = skipped = 0
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(args.vcf, "rt") as source, args.out.open("w") as destination:
        for line in source:
            if line.startswith("#CHROM"):
                samples = line.rstrip().split("\t")[9:]
                expected = [f"mainland_{index}" for index in range(1, 44)] + [
                    f"island_{index}" for index in range(1, 25)
                ]
                if samples != expected:
                    raise ValueError("unexpected registered sample order")
                destination.write("##fileformat=VCFv4.2\n")
                destination.write(f"##contig=<ID=1,length={args.end - args.start}>\n")
                destination.write(
                    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                )
                destination.write(line)
                continue
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            position = int(fields[1]) - 1
            if position < args.start:
                continue
            if position >= args.end:
                break
            if len(fields[3]) != 1 or len(fields[4]) != 1 or "," in fields[4]:
                skipped += 1
                continue
            fields[0] = "1"
            fields[1] = str(position - args.start + 1)
            fields[2] = "."
            fields[5] = "."
            fields[6] = "PASS"
            fields[7] = "."
            fields[8] = "GT"
            destination.write("\t".join(fields) + "\n")
            records += 1
    if len(samples) != 67 or records == 0:
        raise ValueError(f"invalid prepared window: samples={len(samples)}, records={records}")
    report = {
        "source": str(args.vcf),
        "output": str(args.out),
        "source_start": args.start,
        "source_end": args.end,
        "sequence_length": args.end - args.start,
        "diploid_samples": len(samples),
        "records": records,
        "skipped_non_biallelic": skipped,
    }
    args.out.with_suffix(args.out.suffix + ".json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
