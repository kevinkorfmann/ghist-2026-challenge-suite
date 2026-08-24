from __future__ import annotations

import argparse
import gzip
from collections import Counter
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(args.input, "rt") as source, gzip.open(args.output, "wt") as destination:
        for line in source:
            if line.startswith("#"):
                destination.write(line)
                continue
            fields = line.rstrip().split("\t")
            if "," in fields[4]:
                continue
            genotypes = [sample.split(":", 1)[0] for sample in fields[9:]]
            counts: Counter[str] = Counter(
                allele for genotype in genotypes for allele in genotype.replace("/", "|").split("|") if allele != "."
            )
            major = max(counts, key=lambda allele: (counts[allele], -int(allele)))
            for index in range(9, len(fields)):
                parts = fields[index].split(":")
                separator = "|" if "|" in parts[0] else "/"
                alleles = [major if allele == "." else allele for allele in parts[0].split(separator)]
                parts[0] = "|".join(alleles)
                fields[index] = ":".join(parts)
            destination.write("\t".join(fields) + "\n")


if __name__ == "__main__":
    main()
