from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    samples: list[str] = []
    contigs: dict[str, int | None] = {}
    variants = 0
    genotypes = 0
    missing = 0
    phased = 0
    unphased = 0
    multiallelic = 0
    per_contig: Counter[str] = Counter()
    with gzip.open(args.vcf, "rt") as handle:
        for line in handle:
            if line.startswith("##contig=<"):
                content = line.strip()[10:-1]
                fields = dict(item.split("=", 1) for item in content.split(",") if "=" in item)
                contigs[fields["ID"]] = int(fields["length"]) if "length" in fields else None
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip().split("\t")[9:]
                continue
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            variants += 1
            per_contig[fields[0]] += 1
            multiallelic += "," in fields[4]
            for sample in fields[9:]:
                genotype = sample.split(":", 1)[0]
                genotypes += 1
                if "." in genotype:
                    missing += 1
                elif "|" in genotype:
                    phased += 1
                else:
                    unphased += 1

    summary = {
        "state": "COMPLETED",
        "vcf": {"path": str(args.vcf), "sha256": sha256(args.vcf)},
        "samples": len(samples),
        "sample_ids": samples,
        "contigs": contigs,
        "variants": variants,
        "variants_per_contig": dict(per_contig),
        "multiallelic_records": multiallelic,
        "genotypes": genotypes,
        "missing_genotypes": missing,
        "missing_fraction": missing / max(genotypes, 1),
        "phased_nonmissing": phased,
        "unphased_nonmissing": unphased,
        "submission_count": 0,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
