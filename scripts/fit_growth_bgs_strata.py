from __future__ import annotations

import argparse
import json
from pathlib import Path

from ghist_growth_bgs.spectrum import (
    combine_spectra,
    fit_growth_spectrum,
    spectrum_from_counts,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--projection", type=int, default=84)
    parser.add_argument("--restarts", type=int, default=48)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--maxiter", type=int, default=500)
    parser.add_argument("--mutation-rate", type=float, default=3.36e-8)
    parser.add_argument("--seed", type=int, default=75193299)
    args = parser.parse_args()

    audit = json.loads(args.audit.read_text())
    haplotypes = int(audit["expected_haplotypes"])
    distance = audit["distance_strata"]
    joint = audit["joint_strata"]
    distance_bp = audit["physical_bp_by_distance_stratum"]
    joint_bp = audit["physical_bp_by_joint_stratum"]

    def dense(row: dict):
        return spectrum_from_counts(row["unfolded_sfs"], haplotypes)

    candidates: dict[str, tuple[list, int]] = {
        "whole_genome": (
            [dense(row) for row in distance.values()],
            int(audit["sequence_length"]),
        ),
        "coding": ([dense(distance["coding_overlap"])], int(distance_bp["coding_overlap"])),
        "all_noncoding": (
            [dense(row) for key, row in distance.items() if key.startswith("noncoding_")],
            sum(int(value) for key, value in distance_bp.items() if key.startswith("noncoding_")),
        ),
    }
    for quartile in range(1, 5):
        suffix = f"positive_q{quartile}"
        keys = [
            key
            for key in joint
            if key.startswith("noncoding_") and key.endswith(f"|{suffix}")
        ]
        candidates[f"noncoding_{suffix}"] = (
            [dense(joint[key]) for key in keys],
            sum(int(joint_bp[key]) for key in keys),
        )
    far_keys = [
        key
        for key in joint
        if key.startswith("noncoding_10kb_100kb|")
        or key.startswith("noncoding_gt100kb|")
    ]
    candidates["noncoding_gt10kb"] = (
        [dense(joint[key]) for key in far_keys],
        sum(int(joint_bp[key]) for key in far_keys),
    )

    results = {}
    for index, (label, (parts, length)) in enumerate(candidates.items()):
        result = fit_growth_spectrum(
            combine_spectra(parts),
            sequence_length=length,
            mutation_rate=args.mutation_rate,
            projection=args.projection,
            restarts=args.restarts,
            workers=args.workers,
            maxiter=args.maxiter,
            seed=args.seed + 1009 * index,
        )
        result["stratum"] = label
        results[label] = result
        print(json.dumps({label: result}, indent=2), flush=True)

    payload = {
        "source_audit": str(args.audit),
        "projection": args.projection,
        "mutation_rate": args.mutation_rate,
        "results": results,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n")


if __name__ == "__main__":
    main()
