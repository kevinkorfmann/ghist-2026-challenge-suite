from __future__ import annotations

import argparse
import json
from pathlib import Path

import msprime
import numpy as np

from ghist_secondary_contact.demography import SecondaryContactParameters, build_demography
from ghist_secondary_contact.simulate import SimulationNuisance, nuisance_size_changes
from ghist_secondary_contact.spectrum import fit_secondary_contact_spectrum


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--replicate", type=int, required=True)
    parser.add_argument("--truth", type=float, nargs=6, required=True)
    parser.add_argument("--blocks", type=int, default=10)
    parser.add_argument("--block-length", type=float, default=10_000_000)
    parser.add_argument("--projection", type=int, nargs=2, default=(40, 30))
    parser.add_argument("--restarts", type=int, default=12)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--maxiter", type=int, default=400)
    parser.add_argument("--seed", type=int, default=20260919)
    args = parser.parse_args()

    ancestral, size1, size2, split, contact, migration = args.truth
    parameters = SecondaryContactParameters(
        ancestral_size=ancestral,
        population1_size=size1,
        population2_size=size2,
        split_time=split,
        contact_time=contact,
        migration_rate=migration,
    )
    nuisance = SimulationNuisance(
        ancestral_epoch_factors=(0.5, 1.5, 2.0),
        population1_epoch_factors=(0.8, 0.5),
        population2_epoch_factor=1.5,
        population2_contact_factor=1 / 3,
    )
    demography = build_demography(
        parameters,
        size_changes=nuisance_size_changes(parameters, nuisance),
    )
    spectrum = np.zeros((87, 49), dtype=np.float64)
    retained = multiallelic = 0
    rng = np.random.default_rng(args.seed + args.replicate * 1_000_003)
    for _ in range(args.blocks):
        ancestry = msprime.sim_ancestry(
            samples={"population1": 43, "population2": 24},
            ploidy=2,
            demography=demography,
            sequence_length=args.block_length,
            recombination_rate=5.39e-8,
            random_seed=int(rng.integers(1, 2**31 - 2)),
        )
        mutated = msprime.sim_mutations(
            ancestry,
            rate=8.4e-8,
            model=msprime.JC69(),
            discrete_genome=True,
            random_seed=int(rng.integers(1, 2**31 - 2)),
        )
        for variant in mutated.variants():
            if len(variant.alleles) != 2:
                multiallelic += 1
                continue
            genotypes = variant.genotypes
            count1 = int(genotypes[:86].sum())
            count2 = int(genotypes[86:].sum())
            if count1 in (0, 86) and count2 in (0, 48):
                continue
            spectrum[count1, count2] += 1
            retained += 1

    result = fit_secondary_contact_spectrum(
        spectrum,
        mutation_rate=8.4e-8,
        sequence_length=args.blocks * args.block_length,
        restarts=args.restarts,
        workers=args.workers,
        projection=tuple(args.projection),
        maxiter=args.maxiter,
        seed=args.seed + args.replicate * 7_919,
        folded=False,
        convention="generator_matched",
    )
    truth = np.asarray((size1, size2, split, contact, migration))
    estimate = np.asarray(
        (
            result["population1_effective_size"],
            result["population2_effective_size"],
            result["generations_since_split"],
            result["generations_since_contact"],
            result["migration_rate"],
        )
    )
    result.update(
        {
            "replicate": args.replicate,
            "truth": truth.tolist(),
            "ancestral_truth": ancestral,
            "relative_error": ((estimate - truth) / truth).tolist(),
            "competition_rrmse": float(
                np.sqrt(np.sum(((estimate - truth) / truth) ** 2))
            ),
            "retained_biallelic_sites": retained,
            "filtered_multiallelic_sites": multiallelic,
            "blocks": args.blocks,
            "block_length": args.block_length,
            "mutation_model": "JC69",
            "discrete_genome": True,
        }
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps(
            {
                key: value
                for key, value in result.items()
                if key not in {"failed_restarts"}
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
