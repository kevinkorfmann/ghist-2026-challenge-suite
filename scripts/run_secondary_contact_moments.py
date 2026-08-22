from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from ghist_secondary_contact.spectrum import fit_secondary_contact_spectrum
from ghist_secondary_contact.vcf_fast import scan_secondary_contact_vcf


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--restarts", type=int, default=64)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--projection", type=int, nargs=2)
    parser.add_argument("--maxiter", type=int, default=750)
    parser.add_argument("--seed", type=int, default=74031)
    parser.add_argument("--folded", action="store_true")
    parser.add_argument(
        "--convention", choices=("simple", "generator_matched"), default="simple"
    )
    parser.add_argument("--mutation-rate", type=float, default=8.4e-8)
    parser.add_argument("--sequence-length", type=float, default=100_000_000)
    parser.add_argument("--integration-dt-fac", type=float, default=0.02)
    parser.add_argument("--adaptive-integration", action="store_true")
    parser.add_argument("--initial-scaled", type=float, nargs=5)
    args = parser.parse_args()
    spectrum, scan_metadata = scan_secondary_contact_vcf(args.vcf)
    result = fit_secondary_contact_spectrum(
        spectrum,
        mutation_rate=args.mutation_rate,
        sequence_length=args.sequence_length,
        restarts=args.restarts,
        workers=args.workers,
        projection=None if args.projection is None else tuple(args.projection),
        maxiter=args.maxiter,
        seed=args.seed,
        folded=args.folded,
        convention=args.convention,
        integration_dt_fac=args.integration_dt_fac,
        adaptive_integration=args.adaptive_integration,
        initial_scaled=args.initial_scaled,
    )
    result["sfs"] = scan_metadata
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    np.savez_compressed(
        args.out.with_suffix(".sfs.npz"),
        spectrum=spectrum,
        population_names=np.asarray(["mainland", "island"]),
        projection=np.asarray(spectrum.shape, dtype=np.int64) - 1,
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
