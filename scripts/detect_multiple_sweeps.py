from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from ghist_single_sweep.multiple import calibrate_maximum, detect_regions
from ghist_single_sweep.peaks import read_raisd
from scripts.select_background_sweep_candidates import read_exons


def fingerprint(path):
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def main():
    parser = argparse.ArgumentParser(
        description="Generate variable-count candidates from independently calibrated sweep scans"
    )
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument(
        "--null-report",
        type=Path,
        action="append",
        required=True,
        help="Independent matched whole-chromosome null scan; repeat for each replicate",
    )
    parser.add_argument(
        "--context",
        type=Path,
        required=True,
        help="JSON describing matching simulation and scan settings; copied into provenance",
    )
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--exons", type=Path, required=True)
    parser.add_argument("--chromosome", required=True)
    parser.add_argument("--chromosome-length", type=int, required=True)
    parser.add_argument("--separation", type=int, required=True)
    parser.add_argument("--half-width", type=int, required=True)
    parser.add_argument("--max-exon-snap", type=int, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    context = json.loads(args.context.read_text())
    if (
        context.get("demography") != "bottleneck_growth"
        or context.get("background_selection") is not True
    ):
        raise ValueError(
            "Challenge 8 requires matched bottleneck/growth and background-selection nulls"
        )
    paths = [args.report, *args.null_report]
    fingerprints = [fingerprint(path) for path in paths]
    if len({item["sha256"] for item in fingerprints}) != len(paths):
        raise ValueError(
            "Observed and null reports must be distinct; duplicate nulls are not independent replicates"
        )
    maxima = [float(read_raisd(path)[:, 6].max()) for path in args.null_report]
    threshold = calibrate_maximum(maxima, args.alpha)
    result = detect_regions(
        read_raisd(args.report),
        threshold=threshold,
        separation=args.separation,
        half_width=args.half_width,
        exons=read_exons(args.exons, args.chromosome),
        chromosome_length=args.chromosome_length,
        max_exon_snap=args.max_exon_snap,
    )
    result.update(
        {
            "status": "candidate; matched recovery required before freeze",
            "chromosome": args.chromosome.removeprefix("chr"),
            "calibration": {
                "threshold": threshold,
                "alpha": args.alpha,
                "null_maxima": maxima,
                "context": context,
                "assumption": "Observed null and simulated chromosome maxima are exchangeable",
            },
            "localization": {
                key: getattr(args, key)
                for key in (
                    "separation",
                    "half_width",
                    "max_exon_snap",
                    "chromosome_length",
                )
            },
            "source_artifacts": fingerprints
            + [fingerprint(args.exons), fingerprint(args.context)],
            "public_testing_answers_used": False,
        }
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps(
            {key: result[key] for key in ("status", "peak_count", "region_count")},
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
