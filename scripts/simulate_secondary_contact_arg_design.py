from __future__ import annotations

import argparse
import json
from pathlib import Path

import msprime
import numpy as np

from ghist_secondary_contact.demography import SecondaryContactParameters, build_demography
from ghist_secondary_contact.simulate import (
    draw_nuisance,
    focused_parameter_design,
    nuisance_size_changes,
)
from ghist_secondary_contact.features import summarize_tree_sequence


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prefix", type=Path, required=True)
    parser.add_argument("--cell-id", type=int, required=True)
    parser.add_argument("--points", type=int, default=768)
    parser.add_argument("--replicates", type=int, default=3)
    parser.add_argument("--center", type=float, nargs=6, required=True)
    parser.add_argument("--sequence-length", type=float, default=100_000)
    parser.add_argument("--mutation-rate", type=float, default=8.4e-8)
    parser.add_argument("--recombination-rate", type=float, default=5.39e-8)
    parser.add_argument("--seed", type=int, default=20260903)
    args = parser.parse_args()
    total = args.points * args.replicates
    if not 0 <= args.cell_id < total:
        raise ValueError(f"cell id must lie within [0, {total})")
    points, partitions, domains = focused_parameter_design(
        args.points,
        np.asarray(args.center, dtype=np.float64),
        seed=args.seed,
        include_full_prior=False,
    )
    point_id, replicate = divmod(args.cell_id, args.replicates)
    ancestral, mainland, island, split, contact, migration = points[point_id]
    parameters = SecondaryContactParameters(
        ancestral_size=ancestral,
        population1_size=mainland,
        population2_size=island,
        split_time=split,
        contact_time=contact,
        migration_rate=migration,
    )
    rng = np.random.default_rng(args.seed + 1_000_003 * args.cell_id)
    nuisance = draw_nuisance(rng)
    ancestry_seed = int(rng.integers(1, 2**31 - 2))
    mutation_seed = int(rng.integers(1, 2**31 - 2))
    ancestry = msprime.sim_ancestry(
        samples={"population1": 43, "population2": 24},
        ploidy=2,
        demography=build_demography(
            parameters,
            size_changes=nuisance_size_changes(parameters, nuisance),
            migration_asymmetry=nuisance.migration_asymmetry,
        ),
        sequence_length=args.sequence_length,
        recombination_rate=args.recombination_rate * nuisance.recombination_multiplier,
        random_seed=ancestry_seed,
    )
    mutated = msprime.sim_mutations(
        ancestry,
        rate=args.mutation_rate * nuisance.mutation_multiplier,
        model=msprime.BinaryMutationModel(),
        discrete_genome=False,
        random_seed=mutation_seed,
    )
    bundle = summarize_tree_sequence(
        mutated,
        mutated.samples(population=1),
        mutated.samples(population=2),
        n_windows=1,
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
    nuisance_values = np.asarray(
        (
            nuisance.mutation_multiplier,
            nuisance.recombination_multiplier,
            nuisance.migration_asymmetry,
            *(nuisance.ancestral_epoch_factors or (np.nan,) * 3),
            *(nuisance.population1_epoch_factors or (np.nan,) * 2),
            np.nan
            if nuisance.population2_epoch_factor is None
            else nuisance.population2_epoch_factor,
        ),
        dtype=np.float64,
    )
    np.savez_compressed(
        args.prefix.with_suffix(".features.npz"),
        cell_id=args.cell_id,
        point_id=point_id + 2_000_000,
        replicate=replicate,
        partition=partitions[point_id],
        design_domain=domains[point_id],
        truth=np.asarray((mainland, island, split, contact, migration)),
        ancestral_size=ancestral,
        nuisance=nuisance_values,
        truth_hazards=bundle.global_summary[-96:].reshape(3, 32),
        n_variants=bundle.n_variants,
    )
    report = {
        "cell_id": args.cell_id,
        "point_id": point_id + 2_000_000,
        "replicate": replicate,
        "partition": int(partitions[point_id]),
        "design_domain": int(domains[point_id]),
        "parameters": {
            "ancestral_size": ancestral,
            "mainland_size": mainland,
            "island_size": island,
            "split_time": split,
            "contact_time": contact,
            "migration_rate": migration,
        },
        "sites": mutated.num_sites,
        "trees": mutated.num_trees,
        "ancestry_seed": ancestry_seed,
        "mutation_seed": mutation_seed,
    }
    args.prefix.with_suffix(".simulation.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
