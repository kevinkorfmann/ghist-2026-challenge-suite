from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from ghist_secondary_contact.spectrum import fit_secondary_contact_spectrum


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--convention", choices=("simple", "generator_matched"))
    parser.add_argument("--projection", type=int, nargs=2)
    parser.add_argument("--restarts", type=int, default=24)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--maxiter", type=int, default=500)
    parser.add_argument("--seed", type=int, default=20260911)
    parser.add_argument("--folded", action="store_true")
    args = parser.parse_args()
    archive = np.load(args.archive)
    result = fit_secondary_contact_spectrum(
        archive["spectrum"],
        mutation_rate=8.4e-8,
        sequence_length=float(archive["sequence_length"]),
        restarts=args.restarts,
        workers=args.workers,
        projection=None if args.projection is None else tuple(args.projection),
        maxiter=args.maxiter,
        seed=args.seed,
        folded=args.folded,
        convention=args.convention,
    )
    truth = archive["truth"]
    estimate = np.asarray(
        (
            result["population1_effective_size"],
            result["population2_effective_size"],
            result["generations_since_split"],
            result["generations_since_contact"],
            result["migration_rate"],
        )
    )
    result["truth"] = truth.tolist()
    result["relative_error"] = ((estimate - truth) / truth).tolist()
    result["competition_rrmse"] = float(np.sqrt(np.sum(((estimate - truth) / truth) ** 2)))
    result["ancestral_truth"] = float(archive["ancestral_size"])
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({
        key: value for key, value in result.items()
        if key not in {"sfs", "failed_restarts"}
    }, indent=2))


if __name__ == "__main__":
    main()
