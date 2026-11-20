from __future__ import annotations

import argparse
import json
from pathlib import Path

import moments
import msprime
import numpy as np

from ghist_secondary_contact.demography import (
    SecondaryContactParameters,
    build_demography,
)
from ghist_secondary_contact.simulate import SimulationNuisance, nuisance_size_changes
from ghist_secondary_contact.spectrum import (
    generator_matched_secondary_contact_spectrum,
    secondary_contact_spectrum,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--replicates", type=int, default=20)
    parser.add_argument("--sequence-length", type=float, default=500_000)
    parser.add_argument("--seed", type=int, default=20260907)
    args = parser.parse_args()
    parameters = SecondaryContactParameters(
        ancestral_size=20_000,
        population1_size=40_000,
        population2_size=12_000,
        split_time=8_000,
        contact_time=1_000,
        migration_rate=1e-4,
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
    sites = 0
    for replicate in range(args.replicates):
        ancestry = msprime.sim_ancestry(
            samples={"population1": 43, "population2": 24},
            ploidy=2,
            demography=demography,
            sequence_length=args.sequence_length,
            recombination_rate=5.39e-8,
            random_seed=args.seed + 2 * replicate,
        )
        mutated = msprime.sim_mutations(
            ancestry,
            rate=8.4e-8,
            model=msprime.BinaryMutationModel(),
            discrete_genome=False,
            random_seed=args.seed + 2 * replicate + 1,
        )
        for variant in mutated.variants():
            genotypes = variant.genotypes
            count1 = int(genotypes[:86].sum())
            count2 = int(genotypes[86:].sum())
            if count1 + count2 in (0, 134):
                continue
            spectrum[count1, count2] += 1
            sites += 1
    data = moments.Spectrum(spectrum)
    scaled_migration = 2 * parameters.ancestral_size * parameters.migration_rate
    generator_expected = generator_matched_secondary_contact_spectrum(
        (
            parameters.population1_size / parameters.ancestral_size,
            parameters.population2_size / parameters.ancestral_size,
            parameters.split_time / (2 * parameters.ancestral_size),
            parameters.contact_time / parameters.split_time,
            scaled_migration,
        ),
        data.sample_sizes,
    )
    simple_expected = secondary_contact_spectrum(
        (
            parameters.population1_size / parameters.ancestral_size,
            parameters.population2_size / parameters.ancestral_size,
            (parameters.split_time - parameters.contact_time)
            / (2 * parameters.ancestral_size),
            parameters.contact_time / (2 * parameters.ancestral_size),
            scaled_migration,
        ),
        data.sample_sizes,
    )
    generator_likelihood = float(
        moments.Inference.ll_multinom(generator_expected, data)
    )
    simple_likelihood = float(moments.Inference.ll_multinom(simple_expected, data))
    result = {
        "replicates": args.replicates,
        "sequence_length_per_replicate": args.sequence_length,
        "segregating_sites": sites,
        "generator_matched_log_likelihood": generator_likelihood,
        "simple_log_likelihood": simple_likelihood,
        "generator_matched_log_likelihood_gain": generator_likelihood
        - simple_likelihood,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    np.savez_compressed(
        args.out.with_suffix(".sfs.npz"),
        spectrum=spectrum,
        truth=np.asarray(
            (
                parameters.population1_size,
                parameters.population2_size,
                parameters.split_time,
                parameters.contact_time,
                parameters.migration_rate,
            )
        ),
        ancestral_size=parameters.ancestral_size,
        sequence_length=args.replicates * args.sequence_length,
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
