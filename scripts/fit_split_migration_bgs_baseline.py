from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from ghist_split_migration_bgs.spectrum import (
    fit_isolation_with_migration_spectrum,
    stratified_joint_spectra_from_vcf,
)


DISTANCE_EXPOSURE_BP = {
    "coding_overlap": 18_664_670,
    "noncoding_1bp_1kb": 25_925_351,
    "noncoding_1kb_10kb": 41_224_691,
    "noncoding_10kb_100kb": 20_740_375,
    "noncoding_gt100kb": 2_435_119,
}

STRATA = {
    "coding": ("coding_overlap",),
    "near": ("noncoding_1bp_1kb",),
    "middle": ("noncoding_1kb_10kb", "noncoding_10kb_100kb"),
    "far": ("noncoding_gt100kb",),
    "all_noncoding": (
        "noncoding_1bp_1kb",
        "noncoding_1kb_10kb",
        "noncoding_10kb_100kb",
        "noncoding_gt100kb",
    ),
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--coding-bed", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--projection", nargs=2, type=int, default=(60, 50))
    parser.add_argument("--stratum", choices=sorted(STRATA), default="middle")
    parser.add_argument("--restarts", type=int, default=64)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--integration-dt-fac", type=float, default=0.02)
    parser.add_argument("--adaptive-integration", action="store_true")
    parser.add_argument("--initial-fit", type=Path)
    args = parser.parse_args()

    spectra, metadata = stratified_joint_spectra_from_vcf(args.vcf, args.coding_bed)
    included = STRATA[args.stratum]
    selected = np.sum([spectra[label] for label in included], axis=0)
    exposure = sum(DISTANCE_EXPOSURE_BP[label] for label in included)
    initial_scaled = None
    if args.initial_fit is not None:
        scaled = json.loads(args.initial_fit.read_text())["fit"]["scaled_parameters"]
        initial_scaled = (
            scaled["AFR_size_ratio"],
            scaled["EUR_size_ratio"],
            scaled["split_duration_2N"],
            scaled["migration_2N"],
        )
    fit = fit_isolation_with_migration_spectrum(
        selected,
        mutation_rate=1.015e-8,
        sequence_length=exposure,
        projection=tuple(args.projection),
        restarts=args.restarts,
        workers=args.workers,
        integration_dt_fac=args.integration_dt_fac,
        adaptive_integration=args.adaptive_integration,
        initial_scaled=initial_scaled,
    )
    payload = {
        "state": "COMPLETED",
        "stratum": args.stratum,
        "included_distance_labels": list(included),
        "physical_exposure_bp": exposure,
        "metadata": metadata,
        "fit": fit,
        "submission_count": 0,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
