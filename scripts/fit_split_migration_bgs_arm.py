from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from ghist_split_migration_bgs.spectrum import fit_isolation_with_migration_spectrum


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--blocks", type=Path, required=True)
    parser.add_argument("--initial-fit", type=Path, required=True)
    parser.add_argument("--arm", required=True)
    parser.add_argument("--mode", choices=("arm_only", "leave_one_out"), required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--restarts", type=int, default=32)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    archive = np.load(args.blocks)
    contigs = archive["contig"]
    selected = contigs == args.arm
    if args.mode == "leave_one_out":
        selected = ~selected
    if not np.any(selected):
        raise ValueError("arm selection retained no blocks")
    initial = json.loads(args.initial_fit.read_text())["fit"]["scaled_parameters"]
    fit = fit_isolation_with_migration_spectrum(
        np.sum(archive["spectra"][selected], axis=0),
        mutation_rate=1.015e-8,
        sequence_length=float(np.sum(archive["exposure_bp"][selected])),
        projection=(50, 40),
        restarts=args.restarts,
        workers=args.workers,
        maxiter=700,
        seed=75200000 + list(("2L", "2R", "3L", "3R")).index(args.arm),
        initial_scaled=(
            initial["AFR_size_ratio"],
            initial["EUR_size_ratio"],
            initial["split_duration_2N"],
            initial["migration_2N"],
        ),
    )
    payload = {
        "state": "COMPLETED",
        "arm": args.arm,
        "mode": args.mode,
        "blocks": int(np.sum(selected)),
        "physical_exposure_bp": int(np.sum(archive["exposure_bp"][selected])),
        "fit": fit,
        "submission_count": 0,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()

