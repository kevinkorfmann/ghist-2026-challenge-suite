from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from ghist_split_migration_bgs.spectrum import fit_isolation_with_migration_spectrum


PARAMETERS = (
    "AFR_population_size",
    "EUR_population_size",
    "generations_since_split",
    "migration_rate",
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--restarts", type=int, default=64)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()

    files = sorted(
        args.input_dir.glob("task_*.json"),
        key=lambda path: int(path.stem.rsplit("_", 1)[1]),
    )
    if len(files) != 4:
        raise ValueError(f"expected four chromosome-arm simulations, found {len(files)}")
    rows = [json.loads(path.read_text()) for path in files]
    spectra = [np.load(path.with_suffix(".sfs.npz"))["spectrum"] for path in files]
    truth = rows[0]["truth_in_analysis_units"]
    if any(row["truth_in_analysis_units"] != truth for row in rows[1:]):
        raise ValueError("truth differs among chromosome-arm simulations")
    exposure = sum(row["spectrum"]["physical_exposure_bp"] for row in rows)
    fit = fit_isolation_with_migration_spectrum(
        np.sum(spectra, axis=0),
        mutation_rate=1.015e-8,
        sequence_length=exposure,
        projection=(60, 50),
        restarts=args.restarts,
        workers=args.workers,
        maxiter=750,
        seed=75194999,
        initial_scaled=(
            truth["AFR_population_size"] / truth["ancestral_population_size"],
            truth["EUR_population_size"] / truth["ancestral_population_size"],
            truth["generations_since_split"]
            / (2 * truth["ancestral_population_size"]),
            2 * truth["ancestral_population_size"] * truth["migration_rate"],
        ),
    )
    payload = {
        "truth": truth,
        "fit": fit,
        "relative_error": {key: fit[key] / truth[key] - 1 for key in PARAMETERS},
        "chromosome_arms": sorted(row["generator"]["chromosome"] for row in rows),
        "physical_exposure_bp": exposure,
        "accepted_sites": sum(row["spectrum"]["accepted_sites"] for row in rows),
        "submission_count": 0,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
