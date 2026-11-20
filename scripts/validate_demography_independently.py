"""Independent msprime-to-moments recovery diagnostics; never reads challenge data."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import json
import math
from pathlib import Path
import time

import msprime

from ghist_growth_bgs.spectrum import fit_growth_spectrum
from ghist_secondary_contact.demography import (
    SecondaryContactParameters,
    SizeChange,
    build_demography,
)
from ghist_secondary_contact.spectrum import fit_secondary_contact_spectrum
from ghist_split_migration_bgs.spectrum import fit_isolation_with_migration_spectrum


def configurations():
    """Prespecified independent points, not estimated from observed challenge data."""
    rows = []
    for point, parameters in enumerate(
        (
            (12000, 18000, 7000, 10500, 2700, 6e-5),
            (9000, 6500, 21000, 18000, 6500, 1.2e-4),
        )
    ):
        for condition in (
            ("matched", "hidden_epochs", "asymmetric") if point == 0 else ("matched",)
        ):
            rows.append(
                {
                    "task": 2,
                    "point": point,
                    "condition": condition,
                    "parameters": parameters,
                }
            )
    for point, parameters in enumerate(((7000, 9200, 3.7), (15000, 3100, 1.6))):
        for condition in ("exponential", "sudden"):
            rows.append(
                {
                    "task": 3,
                    "point": point,
                    "condition": condition,
                    "parameters": parameters,
                }
            )
    for point, parameters in enumerate(
        ((11000, 8000, 19000, 7800, 7e-5), (18000, 27000, 12000, 12500, 2e-4))
    ):
        for condition in ("matched", "asymmetric") if point == 0 else ("matched",):
            rows.append(
                {
                    "task": 4,
                    "point": point,
                    "condition": condition,
                    "parameters": parameters,
                }
            )
    return rows


def run_one(item):
    config, replicate, out, length, restarts = item
    task, condition = config["task"], config["condition"]
    index = configurations().index(config)
    seed = 202611210 + 1000 * index + 10 * replicate
    started = time.monotonic()
    if task == 2:
        p = SecondaryContactParameters(*config["parameters"])
        changes = (
            []
            if condition != "hidden_epochs"
            else [
                SizeChange(
                    "population1", 0.6 * p.split_time, 0.45 * p.population1_size
                ),
                SizeChange(
                    "population2", 0.4 * p.contact_time, 1.8 * p.population2_size
                ),
                SizeChange("ancestral", 1.4 * p.split_time, 0.65 * p.ancestral_size),
            ]
        )
        demography = build_demography(
            p,
            size_changes=changes,
            migration_asymmetry=4 if condition == "asymmetric" else 1,
        )
        samples = {"population1": 12, "population2": 10}
        mutation_rate, recombination_rate = 8.4e-8, 5.39e-8
        truth = dict(
            zip(
                (
                    "ancestral_effective_size",
                    "population1_effective_size",
                    "population2_effective_size",
                    "generations_since_split",
                    "generations_since_contact",
                    "migration_rate",
                ),
                config["parameters"],
            )
        )
    elif task == 3:
        ancestral, generations, ratio = config["parameters"]
        demography = msprime.Demography()
        demography.add_population(
            name="population",
            initial_size=ancestral * ratio,
            growth_rate=math.log(ratio) / generations
            if condition == "exponential"
            else 0,
        )
        demography.add_population_parameters_change(
            time=generations,
            initial_size=ancestral,
            growth_rate=0,
            population="population",
        )
        samples = {"population": 20}
        mutation_rate, recombination_rate = 3.36e-8, 2e-8
        truth = {
            "ancestral_effective_size": ancestral,
            "generations": generations,
            "post_growth_scaling": ratio,
        }
    else:
        ancestral, first, second, generations, migration = config["parameters"]
        demography = msprime.Demography()
        for name, size in zip(("ANC", "AFR", "EUR"), (ancestral, first, second)):
            demography.add_population(name=name, initial_size=size)
        ratio = 2 if condition == "asymmetric" else 1
        demography.set_migration_rate("AFR", "EUR", migration * ratio)
        demography.set_migration_rate("EUR", "AFR", migration / ratio)
        demography.add_population_split(
            time=generations, derived=["AFR", "EUR"], ancestral="ANC"
        )
        samples = {"AFR": 12, "EUR": 10}
        mutation_rate, recombination_rate = 1.015e-8, 1e-8
        truth = dict(
            zip(
                (
                    "ancestral_population_size",
                    "AFR_population_size",
                    "EUR_population_size",
                    "generations_since_split",
                    "migration_rate",
                ),
                config["parameters"],
            )
        )
    ancestry = msprime.sim_ancestry(
        samples=samples,
        demography=demography,
        sequence_length=length,
        recombination_rate=recombination_rate,
        random_seed=seed,
    )
    mutated = msprime.sim_mutations(
        ancestry,
        rate=mutation_rate,
        random_seed=seed + 1,
        discrete_genome=False,
        model=msprime.BinaryMutationModel(),
    )
    sets = [mutated.samples(population=demography[name].id) for name in samples]
    spectrum = mutated.allele_frequency_spectrum(
        sample_sets=sets, polarised=True, span_normalise=False
    )
    common = {
        "mutation_rate": mutation_rate,
        "sequence_length": length,
        "restarts": restarts,
        "workers": 1,
        "maxiter": 400,
        "seed": seed + 2,
    }
    if task == 2:
        fits = {
            name: fit_secondary_contact_spectrum(
                spectrum,
                projection=(20, 16),
                convention=name,
                initial_scaled=(1.2, 0.9, 0.3, 0.15, 1.0),
                **common,
            )
            for name in ("simple", "generator_matched")
        }
    elif task == 3:
        fits = {
            name: fit_growth_spectrum(
                spectrum, projection=30, convention=name, **common
            )
            for name in ("exponential", "sudden")
        }
    else:
        fits = {
            str(projection): fit_isolation_with_migration_spectrum(
                spectrum,
                projection=projection,
                initial_scaled=(1.2, 0.9, 0.3, 1.0),
                **common,
            )
            for projection in ((20, 16), (12, 10))
        }
    row = {
        "configuration": config,
        "replicate": replicate,
        "seed": seed,
        "truth": truth,
        "fits": fits,
        "relative_error": {
            name: {key: fit[key] / value - 1 for key, value in truth.items()}
            for name, fit in fits.items()
        },
        "sites": mutated.num_sites,
        "length": length,
        "mutation_rate": mutation_rate,
        "recombination_rate": recombination_rate,
        "elapsed_seconds": time.monotonic() - started,
        "public_testing_answers_used": False,
        "initialization": "Fixed generic scaled seed shared across all truths/replicates (Task 2: 1.2,0.9,0.3,0.15,1; Task 4: 1.2,0.9,0.3,1), independent of simulation truth; growth uses global random restarts",
        "scope": "Independent msprime genealogies and infinite-sites mutations; neutral numerical/misspecification diagnostic, not BGS or whole-chromosome coverage validation",
    }
    path = (
        Path(out) / f"task{task}_point{config['point']}_{condition}_rep{replicate}.json"
    )
    path.write_text(json.dumps(row, indent=2) + "\n")
    return str(path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--replicates", type=int, default=2)
    parser.add_argument("--sequence-length", type=int, default=2_000_000)
    parser.add_argument("--restarts", type=int, default=12)
    args = parser.parse_args()
    if min(args.workers, args.replicates, args.sequence_length, args.restarts) < 1:
        raise ValueError("positive execution settings required")
    args.out.mkdir(parents=True, exist_ok=True)
    work = [
        (row, replicate, str(args.out), args.sequence_length, args.restarts)
        for row in configurations()
        for replicate in range(args.replicates)
    ]
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        for path in executor.map(run_one, work):
            print(path, flush=True)


if __name__ == "__main__":
    main()
