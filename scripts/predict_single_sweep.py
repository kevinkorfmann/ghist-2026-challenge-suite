"""Generate and validate a provisional Task 5 candidate without truth input."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from ghist_single_sweep.predict import predict_single_sweep
from ghist_single_sweep.scan import read_bed_intervals, read_phased_biallelic_vcf
from ghist_submission.validate import validate_submission


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def run(vcf: Path, accessible: Path, dataset: str, bed: Path, out: Path) -> dict:
    if any(path.exists() for path in (bed, out)):
        raise ValueError("Refusing to overwrite existing candidate artifacts")
    variants = read_phased_biallelic_vcf(vcf)
    expected = "7" if dataset == "testing" else "9"
    if variants.contig.removeprefix("chr") != expected:
        raise ValueError(f"{dataset} requires chromosome {expected}")
    intervals = read_bed_intervals(accessible, variants.contig)
    if intervals[-1, 1] > variants.contig_length:
        raise ValueError("Accessibility extends beyond the VCF contig")
    result = predict_single_sweep(variants, intervals)
    bed.parent.mkdir(parents=True, exist_ok=True)
    start, end = result["interval"]
    bed.write_text(f"{expected}\t{start}\t{end}\n")
    result["format_validation"] = validate_submission(bed, 5, dataset)
    result["input_provenance"] = {
        "vcf_sha256": sha256(vcf),
        "accessibility_sha256": sha256(accessible),
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2) + "\n")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--accessible", type=Path, required=True)
    parser.add_argument("--dataset", choices=("testing", "final"), required=True)
    parser.add_argument("--bed", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.vcf, args.accessible, args.dataset, args.bed, args.out)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
