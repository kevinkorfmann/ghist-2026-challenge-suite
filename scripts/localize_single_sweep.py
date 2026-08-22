from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from ghist_single_sweep.scan import flank_ld_scan, read_phased_biallelic_vcf


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--region", type=int, nargs=2, required=True)
    parser.add_argument("--grid-step", type=int, default=1000)
    parser.add_argument("--flank-width", type=int, action="append", default=[])
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args()

    flank_widths = args.flank_width or [50_000, 100_000, 250_000, 500_000]
    variants = read_phased_biallelic_vcf(args.vcf)
    candidates = np.arange(args.region[0], args.region[1] + 1, args.grid_step, dtype=np.int64)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    scans: dict[str, dict] = {}
    for flank_width in flank_widths:
        scan = flank_ld_scan(variants, candidates, flank_width)
        path = args.out_dir / f"flank_ld_{flank_width}.csv"
        names = list(scan)
        with path.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(names)
            writer.writerows(zip(*(scan[name] for name in names), strict=True))
        valid = np.isfinite(scan["flank_ld_ratio"])
        order = np.flatnonzero(valid)[np.argsort(scan["flank_ld_ratio"][valid])[::-1]]
        scans[str(flank_width)] = {
            "csv": str(path),
            "maximum": {
                key: float(scan[key][order[0]])
                for key in (
                    "candidate_position",
                    "within_left_r2",
                    "within_right_r2",
                    "between_r2",
                    "flank_ld_ratio",
                )
            },
            "top_positions": [int(scan["candidate_position"][index]) for index in order[:20]],
        }

    left = int(np.searchsorted(variants.positions, args.region[0], side="left"))
    right = int(np.searchsorted(variants.positions, args.region[1], side="right"))
    local_positions = variants.positions[left:right]
    gaps = np.diff(local_positions.astype(np.int64))
    gap_index = int(np.argmax(gaps))
    summary = {
        "state": "COMPLETED",
        "region": args.region,
        "grid_step": args.grid_step,
        "largest_variant_gap": {
            "left": int(local_positions[gap_index]),
            "right": int(local_positions[gap_index + 1]),
            "length": int(gaps[gap_index]),
            "midpoint": float(
                (int(local_positions[gap_index]) + int(local_positions[gap_index + 1])) / 2
            ),
        },
        "flank_ld_scans": scans,
        "submission_count": 0,
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

