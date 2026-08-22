from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pysam

from ghist_secondary_contact.spectrum import (
    fit_secondary_contact_spectrum,
    joint_spectrum_from_vcf,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--restarts", type=int, default=64)
    parser.add_argument("--seed", type=int, default=74031)
    parser.add_argument("--folded", action="store_true")
    parser.add_argument("--mutation-rate", type=float, default=8.4e-8)
    parser.add_argument("--sequence-length", type=float, default=100_000_000)
    args = parser.parse_args()
    vcf = pysam.VariantFile(str(args.vcf))
    samples = list(vcf.header.samples)
    populations = {
        "mainland": [sample for sample in samples if sample.startswith("mainland_")],
        "island": [sample for sample in samples if sample.startswith("island_")],
    }
    spectrum, metadata = joint_spectrum_from_vcf(
        args.vcf,
        populations,
        polarized=not args.folded,
    )
    result = fit_secondary_contact_spectrum(
        spectrum,
        mutation_rate=args.mutation_rate,
        sequence_length=args.sequence_length,
        restarts=args.restarts,
        seed=args.seed,
        folded=args.folded,
    )
    result["sfs"] = metadata
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    np.savez_compressed(
        args.out.with_suffix(".sfs.npz"),
        spectrum=spectrum,
        population_names=np.asarray(metadata["population_names"]),
        projection=np.asarray(metadata["projection"]),
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
