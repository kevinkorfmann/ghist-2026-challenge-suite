from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from ghist_growth_bgs.spectrum import fit_growth_spectrum


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--simulation", type=Path, required=True)
    parser.add_argument("--sfs", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--restarts", type=int, default=32)
    args = parser.parse_args()
    simulation = json.loads(args.simulation.read_text())
    spectrum = np.load(args.sfs)["spectrum"]
    sudden = fit_growth_spectrum(
        spectrum,
        sequence_length=simulation["spectrum"]["physical_exposure_bp"],
        mutation_rate=simulation["truth_in_analysis_units"]["mutation_rate"],
        projection=simulation["spectrum"]["haplotypes"],
        restarts=args.restarts,
        workers=1,
        maxiter=500,
        seed=simulation["generator"]["seed"] + 2,
        convention="sudden",
    )
    exponential = simulation["fit"]
    result = {
        "truth": simulation["truth_in_analysis_units"],
        "exponential": exponential,
        "sudden": sudden,
        "log_likelihood_exponential_minus_sudden": (
            exponential["multinomial_log_likelihood"]
            - sudden["multinomial_log_likelihood"]
        ),
        "preferred": (
            "exponential"
            if exponential["multinomial_log_likelihood"]
            > sudden["multinomial_log_likelihood"]
            else "sudden"
        ),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
