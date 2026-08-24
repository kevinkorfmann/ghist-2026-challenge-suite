from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from ghist_single_sweep.scan import read_bed_intervals, read_phased_biallelic_vcf, scan_windows


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--accessible-bed", type=Path, required=True)
    parser.add_argument("--window-size", type=int, action="append", default=[])
    parser.add_argument(
        "--score-field",
        choices=["composite_score", "soft_composite_score"],
        default="composite_score",
    )
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args()

    window_sizes = args.window_size or [25_000, 50_000, 100_000, 250_000, 500_000]
    variants = read_phased_biallelic_vcf(args.vcf)
    accessibility = read_bed_intervals(args.accessible_bed, variants.contig)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    scale_summaries: dict[str, dict] = {}
    for window_size in window_sizes:
        step_size = max(window_size // 5, 1)
        scan = scan_windows(variants, accessibility, window_size, step_size)
        path = args.out_dir / f"scan_{window_size}.csv"
        names = list(scan)
        with path.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(names)
            writer.writerows(zip(*(scan[name] for name in names), strict=True))
        valid = np.isfinite(scan[args.score_field])
        order = np.flatnonzero(valid)[np.argsort(scan[args.score_field][valid])[::-1]]
        scale_summaries[str(window_size)] = {
            "step_size": step_size,
            "windows": int(valid.sum()),
            "top_windows": [
                {
                    "start": int(scan["start"][index]),
                    "end": int(scan["end"][index]),
                    "center": int(scan["center"][index]),
                    "composite_score": float(scan["composite_score"][index]),
                    "soft_composite_score": float(scan["soft_composite_score"][index]),
                    "pi": float(scan["pi"][index]),
                    "tajima_d": float(scan["tajima_d"][index]),
                    "h12": float(scan["h12"][index]),
                    "h2_h1": float(scan["h2_h1"][index]),
                    "mean_r2": float(scan["mean_r2"][index]),
                }
                for index in order[:20]
            ],
            "csv": str(path),
        }
    summary = {
        "state": "COMPLETED",
        "contig": variants.contig,
        "contig_length": variants.contig_length,
        "diploid_samples": variants.haplotypes.shape[1] // 2,
        "haplotypes": variants.haplotypes.shape[1],
        "accepted_biallelic_sites": variants.positions.size,
        "skipped_records": variants.skipped_records,
        "accessibility_intervals": accessibility.tolist(),
        "files": {
            "vcf": {"path": str(args.vcf), "sha256": sha256(args.vcf)},
            "accessible_bed": {
                "path": str(args.accessible_bed),
                "sha256": sha256(args.accessible_bed),
            },
        },
        "scales": scale_summaries,
        "ranking_score": args.score_field,
        "submission_count": 0,
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
