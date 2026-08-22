from __future__ import annotations

import argparse
import json
from pathlib import Path

import msprime
import numpy as np

from ghist_secondary_contact.demography import SecondaryContactParameters, build_demography


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prefix", type=Path, required=True)
    parser.add_argument("--ancestral-size", type=float, default=18_000)
    parser.add_argument("--mainland-size", type=float, default=24_000)
    parser.add_argument("--island-size", type=float, default=9_000)
    parser.add_argument("--split-time", type=float, default=12_000)
    parser.add_argument("--contact-time", type=float, default=1_500)
    parser.add_argument("--migration-rate", type=float, default=8e-5)
    parser.add_argument("--sequence-length", type=float, default=100_000)
    parser.add_argument("--mutation-rate", type=float, default=8.4e-8)
    parser.add_argument("--recombination-rate", type=float, default=5.39e-8)
    parser.add_argument("--seed", type=int, default=71001)
    args = parser.parse_args()
    parameters = SecondaryContactParameters(
        ancestral_size=args.ancestral_size,
        population1_size=args.mainland_size,
        population2_size=args.island_size,
        split_time=args.split_time,
        contact_time=args.contact_time,
        migration_rate=args.migration_rate,
    )
    ancestry = msprime.sim_ancestry(
        samples={"population1": 43, "population2": 24},
        ploidy=2,
        demography=build_demography(parameters),
        sequence_length=args.sequence_length,
        recombination_rate=args.recombination_rate,
        random_seed=args.seed,
    )
    mutated = msprime.sim_mutations(
        ancestry,
        rate=args.mutation_rate,
        model=msprime.BinaryMutationModel(),
        discrete_genome=False,
        random_seed=args.seed + 1,
    )
    args.prefix.parent.mkdir(parents=True, exist_ok=True)
    mutated.dump(args.prefix.with_suffix(".truth.trees"))
    names = [f"mainland_{index}" for index in range(1, 44)] + [
        f"island_{index}" for index in range(1, 25)
    ]
    with args.prefix.with_suffix(".vcf").open("w") as handle:
        mutated.write_vcf(
            handle,
            individual_names=names,
            position_transform=lambda position: np.asarray(position) + 1,
        )
    report = {
        "prefix": str(args.prefix),
        "parameters": vars(args) | {"prefix": str(args.prefix)},
        "samples": mutated.num_samples,
        "sites": mutated.num_sites,
        "mutations": mutated.num_mutations,
        "trees": mutated.num_trees,
    }
    args.prefix.with_suffix(".simulation.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
