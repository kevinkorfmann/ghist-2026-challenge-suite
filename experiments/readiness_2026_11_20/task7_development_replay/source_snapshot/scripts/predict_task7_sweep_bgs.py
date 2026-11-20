"""Run the new masked ten-estimator Task 7 pipeline; never freeze or submit."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
from pathlib import Path

import numpy as np
import pysam

from ghist_single_sweep.background import SETTINGS, masked_variants, predict_masked_ten
from ghist_single_sweep.peaks import read_raisd
from ghist_single_sweep.scan import read_bed_intervals, read_phased_biallelic_vcf
from ghist_submission.validate import validate_submission

ROOT = Path(__file__).resolve().parents[1]
SOURCE_PATHS = (
    "scripts/predict_task7_sweep_bgs.py",
    "src/ghist_single_sweep/background.py",
    "src/ghist_single_sweep/scan.py",
    "src/ghist_single_sweep/peaks.py",
    "src/ghist_submission/validate.py",
)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def source_manifest():
    return {path: sha256(ROOT / path) for path in SOURCE_PATHS}


def write_masked_vcf(variants, path):
    """Canonical genotype-only export; sample order and physical exposure retained."""
    n = variants.haplotypes.shape[1]
    with path.open("w") as handle:
        handle.write(
            f"##fileformat=VCFv4.2\n##contig=<ID={variants.contig},length={variants.contig_length}>\n"
        )
        handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        handle.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
            + "\t".join(f"sample_{i}" for i in range(n // 2))
            + "\n"
        )
        for pos, row in zip(variants.positions, variants.haplotypes, strict=True):
            handle.write(
                f"{variants.contig}\t{int(pos) + 1}\t.\tA\tC\t.\tPASS\t.\tGT\t"
                + "\t".join(f"{row[i]}|{row[i + 1]}" for i in range(0, n, 2))
                + "\n"
            )
    reread = read_phased_biallelic_vcf(path)
    if (
        reread.contig_length != variants.contig_length
        or not np.array_equal(reread.positions, variants.positions)
        or not np.array_equal(reread.haplotypes, variants.haplotypes)
    ):
        raise ValueError("Canonical masked VCF changed genotype or exposure contract")


def validate_raisd_site_count(info_path, expected):
    counts = [
        int(value) for value in re.findall(r"\|\s+snps\s+(\d+)", info_path.read_text())
    ]
    if counts != [expected]:
        raise ValueError("RAiSD and masked predictor usable SNP counts differ")


def bound_provenance(vcf, accessible, raisd):
    return {
        "input_provenance": {
            "vcf_sha256": sha256(vcf),
            "accessibility_sha256": sha256(accessible),
        },
        "source_sha256": source_manifest(),
        "software": {
            "numpy": np.__version__,
            "pysam": pysam.__version__,
            "raisd_binary_sha256": sha256(raisd),
        },
    }


def run(vcf, accessible, dataset, raisd, out_dir):
    if out_dir.exists():
        raise ValueError("Refusing to overwrite an existing Task 7 output directory")
    bound = bound_provenance(vcf, accessible, raisd)
    expected = "7" if dataset == "testing" else "9"
    original = read_phased_biallelic_vcf(vcf)
    if original.contig.removeprefix("chr") != expected:
        raise ValueError(f"{dataset} requires chromosome {expected}")
    variants, mask = masked_variants(
        original, read_bed_intervals(accessible, original.contig)
    )
    if len(variants.positions) < SETTINGS["raisd_snp_window"]:
        raise ValueError("Too few callable SNPs for RAiSD's 50-SNP window")
    out_dir.mkdir(parents=True)
    canonical = out_dir / "masked.vcf"
    write_masked_vcf(variants, canonical)
    mask_path = out_dir / "accessible.bed"
    mask_path.write_text("".join(f"{expected}\t{int(a)}\t{int(b)}\n" for a, b in mask))
    canonical_hash = sha256(canonical)
    with (out_dir / "raisd.log").open("w") as log:
        subprocess.run(
            [
                str(raisd.resolve()),
                "-n",
                "task7",
                "-I",
                str(canonical.resolve()),
                "-w",
                str(SETTINGS["raisd_snp_window"]),
                "-R",
            ],
            cwd=out_dir,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )
    report_path, info_path = (
        out_dir / "RAiSD_Report.task7",
        out_dir / "RAiSD_Info.task7",
    )
    validate_raisd_site_count(info_path, len(variants.positions))
    result = predict_masked_ten(variants, mask, read_raisd(report_path))
    bed = out_dir / "candidate.bed"
    start, end = result["interval"]
    bed.write_text(f"{expected}\t{start}\t{end}\n")
    result["format_validation"] = validate_submission(bed, 7, dataset)
    result["dataset"] = dataset
    if (
        bound_provenance(vcf, accessible, raisd) != bound
        or sha256(canonical) != canonical_hash
    ):
        raise ValueError(
            "Task 7 input, source, canonical VCF or RAiSD binary changed during inference"
        )
    result.update(bound)
    result["artifact_sha256"] = {
        p.name: sha256(p) for p in (canonical, mask_path, report_path, info_path, bed)
    }
    result["raisd_contract"] = (
        "Same mask-filtered A/C genotypes; reject report windows crossing inaccessible bases; native chromosome-wide RAiSD score normalization retained"
    )
    (out_dir / "prediction.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n"
    )
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--accessible", type=Path, required=True)
    parser.add_argument("--dataset", choices=("testing", "final"), required=True)
    parser.add_argument("--raisd", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    print(
        json.dumps(
            run(args.vcf, args.accessible, args.dataset, args.raisd, args.out_dir),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
