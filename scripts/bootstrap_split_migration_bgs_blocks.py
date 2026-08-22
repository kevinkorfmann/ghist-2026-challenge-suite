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
    parser.add_argument("--replicate", type=int, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--restarts", type=int, default=8)
    args = parser.parse_args()
    archive = np.load(args.blocks)
    spectra = archive["spectra"]
    exposure = archive["exposure_bp"]
    contigs = archive["contig"]
    rng = np.random.default_rng(75198000 + args.replicate)
    selected = []
    for contig in np.unique(contigs):
        indices = np.flatnonzero(contigs == contig)
        selected.extend(rng.choice(indices, size=len(indices), replace=True).tolist())
    initial_payload = json.loads(args.initial_fit.read_text())["fit"]
    scaled = initial_payload["scaled_parameters"]
    fit = fit_isolation_with_migration_spectrum(
        np.sum(spectra[selected], axis=0),
        mutation_rate=1.015e-8,
        sequence_length=float(np.sum(exposure[selected])),
        projection=(50, 40),
        restarts=args.restarts,
        workers=1,
        maxiter=600,
        seed=75199000 + args.replicate,
        initial_scaled=(
            scaled["AFR_size_ratio"],
            scaled["EUR_size_ratio"],
            scaled["split_duration_2N"],
            scaled["migration_2N"],
        ),
    )
    payload = {
        "state": "COMPLETED",
        "replicate": args.replicate,
        "selected_blocks": len(selected),
        "physical_exposure_bp": int(np.sum(exposure[selected])),
        "fit": fit,
        "submission_count": 0,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()

