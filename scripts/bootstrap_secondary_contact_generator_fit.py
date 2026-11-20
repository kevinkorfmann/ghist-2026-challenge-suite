from __future__ import annotations

import argparse
import json
from pathlib import Path

import msprime
import numpy as np

from ghist_secondary_contact.demography import (
    SecondaryContactParameters,
    build_demography,
)
from ghist_secondary_contact.simulate import SimulationNuisance, nuisance_size_changes
from ghist_secondary_contact.spectrum import fit_secondary_contact_spectrum
from ghist_secondary_contact.contract import SPECTRUM_CONTRACT


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--replicate", type=int, required=True)
    parser.add_argument("--truth", type=float, nargs=6, required=True)
    parser.add_argument("--blocks", type=int, default=10)
    parser.add_argument("--block-length", type=float, default=10_000_000)
    parser.add_argument("--projection", type=int, nargs=2, default=(40, 30))
    parser.add_argument("--full-projection", action="store_true")
    parser.add_argument("--restarts", type=int, default=12)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--maxiter", type=int, default=400)
    parser.add_argument("--seed", type=int, default=20260919)
    parser.add_argument("--integration-dt-fac", type=float, default=0.02)
    parser.add_argument("--adaptive-integration", action="store_true")
    parser.add_argument("--initial-scaled", type=float, nargs=5)
    parser.add_argument(
        "--recombination-multipliers", type=float, nargs="+", default=(1.0,)
    )
    parser.add_argument(
        "--mutation-rate-multipliers", type=float, nargs="+", default=(1.0,)
    )
    parser.add_argument(
        "--local-size-multipliers", type=float, nargs="+", default=(1.0,)
    )
    parser.add_argument("--migration-asymmetry", type=float, default=1.0)
    parser.add_argument("--polarization-error", type=float, default=0.0)
    args = parser.parse_args()
    multiplier_groups = (
        args.recombination_multipliers,
        args.mutation_rate_multipliers,
        args.local_size_multipliers,
    )
    if any(any(value <= 0 for value in group) for group in multiplier_groups):
        parser.error("rate and local-size multipliers must be positive")
    if args.migration_asymmetry <= 0:
        parser.error("migration asymmetry must be positive")
    if not 0 <= args.polarization_error < 0.5:
        parser.error("polarization error must be in [0, 0.5)")

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
    spectrum = np.zeros((87, 49), dtype=np.float64)
    retained = multiallelic = 0
    rng = np.random.default_rng(args.seed + args.replicate * 1_000_003)
    for block in range(args.blocks):
        local_size = args.local_size_multipliers[
            block % len(args.local_size_multipliers)
        ]
        local_parameters = SecondaryContactParameters(
            ancestral_size=parameters.ancestral_size * local_size,
            population1_size=parameters.population1_size * local_size,
            population2_size=parameters.population2_size * local_size,
            split_time=parameters.split_time,
            contact_time=parameters.contact_time,
            migration_rate=parameters.migration_rate,
        )
        demography = build_demography(
            local_parameters,
            size_changes=nuisance_size_changes(local_parameters, nuisance),
            migration_asymmetry=args.migration_asymmetry,
        )
        recombination_multiplier = args.recombination_multipliers[
            block % len(args.recombination_multipliers)
        ]
        mutation_multiplier = args.mutation_rate_multipliers[
            block % len(args.mutation_rate_multipliers)
        ]
        ancestry = msprime.sim_ancestry(
            samples={"population1": 43, "population2": 24},
            ploidy=2,
            demography=demography,
            sequence_length=args.block_length,
            recombination_rate=5.39e-8 * recombination_multiplier,
            random_seed=int(rng.integers(1, 2**31 - 2)),
        )
        mutated = msprime.sim_mutations(
            ancestry,
            rate=8.4e-8 * mutation_multiplier,
            model=msprime.JC69(),
            discrete_genome=True,
            random_seed=int(rng.integers(1, 2**31 - 2)),
        )
        polarization_rng = np.random.default_rng(
            args.seed + args.replicate * 1_000_003 + block * 65_537 + 17
        )
        for variant in mutated.variants():
            if len(variant.alleles) != 2:
                multiallelic += 1
                continue
            genotypes = variant.genotypes
            count1 = int(genotypes[:86].sum())
            count2 = int(genotypes[86:].sum())
            if (
                args.polarization_error
                and polarization_rng.random() < args.polarization_error
            ):
                count1 = 86 - count1
                count2 = 48 - count2
            if count1 + count2 in (0, 134):
                continue
            spectrum[count1, count2] += 1
            retained += 1

    result = fit_secondary_contact_spectrum(
        spectrum,
        mutation_rate=8.4e-8,
        sequence_length=args.blocks * args.block_length,
        restarts=args.restarts,
        workers=args.workers,
        projection=None if args.full_projection else tuple(args.projection),
        maxiter=args.maxiter,
        seed=args.seed + args.replicate * 7_919,
        folded=False,
        convention="generator_matched",
        integration_dt_fac=args.integration_dt_fac,
        adaptive_integration=args.adaptive_integration,
        initial_scaled=args.initial_scaled,
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
            "spectrum_contract": SPECTRUM_CONTRACT,
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
            "recombination_multipliers": list(args.recombination_multipliers),
            "mutation_rate_multipliers": list(args.mutation_rate_multipliers),
            "local_size_multipliers": list(args.local_size_multipliers),
            "migration_asymmetry": args.migration_asymmetry,
            "polarization_error": args.polarization_error,
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
