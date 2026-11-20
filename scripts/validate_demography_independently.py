"""Independent, resumable msprime-to-moments recovery; no challenge inputs."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
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
from ghist_inference.checkpoint import (
    atomic_json,
    load_spectrum_checkpoint,
    save_spectrum_checkpoint,
)

SIMULATION_PROTOCOL = "independent_msprime_diploid_v2"
FIT_PROTOCOL = "bounded_refinement_diagnostics_v1"


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


def simulation_seeds(configuration_index, replicate, scheme="counter_v2"):
    """Disjoint ancestry, mutation and optimizer streams without hash collisions.

    The fixed radix reserves sixteen configuration IDs and four streams each.
    Increasing the replicate count never changes an existing allocation. The
    old arithmetic is retained only for exact replay of fewer than 101 reps.
    """
    if not isinstance(configuration_index, int) or not 0 <= configuration_index < 16:
        raise ValueError("configuration index must be an integer in [0, 16)")
    if not isinstance(replicate, int) or replicate < 0:
        raise ValueError("replicate must be a nonnegative integer")
    if scheme == "counter_v2":
        seed = 202611210 + 64 * replicate + 4 * configuration_index
        values = (seed, seed + 1, seed + 2, seed + 3)
    elif scheme == "legacy_v1":
        if replicate >= 100:
            raise ValueError(
                "legacy seed allocation collides at replicate 100; use counter_v2"
            )
        seed = 202611210 + 1000 * configuration_index + 10 * replicate
        values = (seed, seed + 1, seed + 2, seed + 2)
    else:
        raise ValueError("unknown seed scheme")
    if max(values) >= 2**32:
        raise ValueError("seed allocation exceeds the msprime uint32 range")
    return dict(zip(("ancestry", "mutations", "fit0", "fit1"), values))


def simulate_configuration(config, seeds, length):
    task, condition = config["task"], config["condition"]
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
        random_seed=seeds["ancestry"],
    )
    mutated = msprime.sim_mutations(
        ancestry,
        rate=mutation_rate,
        random_seed=seeds["mutations"],
        discrete_genome=False,
        model=msprime.BinaryMutationModel(),
    )
    sets = [mutated.samples(population=demography[name].id) for name in samples]
    spectrum = mutated.allele_frequency_spectrum(
        sample_sets=sets, polarised=True, span_normalise=False
    )
    return spectrum, {
        "truth": truth,
        "sites": mutated.num_sites,
        "length": length,
        "mutation_rate": mutation_rate,
        "recombination_rate": recombination_rate,
        "msprime_version": msprime.__version__,
        "coalescent_ploidy": 2,
    }


def run_one(item):
    config, replicate, out, length, restarts, scheme = item
    started = time.monotonic()
    task, condition = config["task"], config["condition"]
    index = configurations().index(config)
    seeds = simulation_seeds(index, replicate, scheme)
    stem = f"task{task}_point{config['point']}_{condition}_rep{replicate}"
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    specification = {
        "simulation_protocol": SIMULATION_PROTOCOL,
        "configuration": config,
        "replicate": replicate,
        "length": length,
        "seeds": seeds,
        "seed_scheme": scheme,
    }
    checkpoint = out / (stem + ".sfs.npz")
    if checkpoint.exists():
        spectrum, metadata = load_spectrum_checkpoint(
            checkpoint, expected=specification
        )
    else:
        spectrum, metadata = simulate_configuration(config, seeds, length)
        metadata = save_spectrum_checkpoint(
            checkpoint, spectrum, dict(metadata, **specification)
        )
    common = {
        "mutation_rate": metadata["mutation_rate"],
        "sequence_length": length,
        "restarts": restarts,
        "workers": 1,
        "maxiter": 400,
    }
    if task == 2:
        models = [
            (
                name,
                {
                    "convention": name,
                    "projection": (20, 16),
                    "initial_scaled": (1.2, 0.9, 0.3, 0.15, 1.0),
                },
            )
            for name in ("simple", "generator_matched")
        ]
        fitter = fit_secondary_contact_spectrum
    elif task == 3:
        models = [
            (name, {"convention": name, "projection": 30})
            for name in ("exponential", "sudden")
        ]
        fitter = fit_growth_spectrum
    else:
        models = [
            (
                str(projection),
                {"projection": projection, "initial_scaled": (1.2, 0.9, 0.3, 1.0)},
            )
            for projection in ((20, 16), (12, 10))
        ]
        fitter = fit_isolation_with_migration_spectrum
    fits, failures = {}, {}
    for model_index, (name, settings) in enumerate(models):
        settings = dict(common, **settings, seed=seeds[f"fit{model_index}"])
        fit_spec = json.loads(
            json.dumps(
                {
                    "fit_protocol": FIT_PROTOCOL,
                    "settings": settings,
                    "spectrum_sha256": metadata["spectrum_sha256"],
                }
            )
        )
        fit_path = out / (stem + f".fit{model_index}.json")
        record = None
        if fit_path.exists():
            record = json.loads(fit_path.read_text())
            if record["specification"] != fit_spec:
                raise ValueError(f"Fit checkpoint configuration mismatch: {fit_path}")
            if record["status"] not in ("completed", "failed"):
                raise ValueError(f"Unknown fit checkpoint status: {fit_path}")
        if record is None or record["status"] == "failed":
            previous_attempts = (
                []
                if record is None
                else [
                    *record.get("previous_failed_attempts", []),
                    {key: record[key] for key in ("error", "elapsed_seconds")},
                ]
            )
            fit_started = time.monotonic()
            try:
                record = {"status": "completed", "fit": fitter(spectrum, **settings)}
            except Exception as error:
                record = {
                    "status": "failed",
                    "error": f"{type(error).__name__}: {error}",
                }
            record.update(
                specification=fit_spec,
                model=name,
                elapsed_seconds=time.monotonic() - fit_started,
                previous_failed_attempts=previous_attempts,
            )
            atomic_json(fit_path, record)
        if record["status"] == "completed":
            fits[name] = record["fit"]
        else:
            failures[name] = record["error"]
    truth = metadata["truth"]
    row = {
        "status": "completed" if not failures else "partial" if fits else "failed",
        "configuration": config,
        "replicate": replicate,
        "seed": seeds["ancestry"],
        "seeds": seeds,
        "seed_scheme": scheme,
        "truth": truth,
        "fits": fits,
        "failures": failures,
        "relative_error": {
            name: {key: fit[key] / value - 1 for key, value in truth.items()}
            for name, fit in fits.items()
        },
        "sites": metadata["sites"],
        "length": length,
        "mutation_rate": metadata["mutation_rate"],
        "recombination_rate": metadata["recombination_rate"],
        "elapsed_seconds_this_invocation": time.monotonic() - started,
        "public_testing_answers_used": False,
        "spectrum_sha256": metadata["spectrum_sha256"],
        "simulation_protocol": SIMULATION_PROTOCOL,
        "fit_protocol": FIT_PROTOCOL,
        "initialization": "Fixed generic scaled starts shared across all truths (Tasks 2/4), independent of simulation truth; growth uses global random starts",
        "scope": f"Independent neutral linked {length} bp chromosomes; not BGS or coverage validation",
    }
    path = out / (stem + ".json")
    atomic_json(path, row)
    return str(path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--replicates", type=int, default=2)
    parser.add_argument("--sequence-length", type=int, default=2_000_000)
    parser.add_argument("--restarts", type=int, default=12)
    parser.add_argument(
        "--seed-scheme", choices=("counter_v2", "legacy_v1"), default="counter_v2"
    )
    parser.add_argument(
        "--tasks", nargs="+", type=int, choices=(2, 3, 4), default=(2, 3, 4)
    )
    args = parser.parse_args()
    if min(args.workers, args.replicates, args.sequence_length, args.restarts) < 1:
        raise ValueError("positive execution settings required")
    selected = [row for row in configurations() if row["task"] in args.tasks]
    work = [
        (
            row,
            replicate,
            str(args.out),
            args.sequence_length,
            args.restarts,
            args.seed_scheme,
        )
        for row in selected
        for replicate in range(args.replicates)
    ]
    # Validate the whole plan before spending compute or creating partial files.
    for row, replicate, *_ in work:
        simulation_seeds(configurations().index(row), replicate, args.seed_scheme)
    args.out.mkdir(parents=True, exist_ok=True)
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        pending = {executor.submit(run_one, item): item for item in work}
        for future in as_completed(pending):
            item = pending[future]
            try:
                print(future.result(), flush=True)
            except Exception as error:
                # A corrupt/incompatible checkpoint must be visible and never overwritten.
                print(
                    json.dumps(
                        {
                            "status": "execution_failed",
                            "configuration": item[0],
                            "replicate": item[1],
                            "error": f"{type(error).__name__}: {error}",
                        }
                    ),
                    flush=True,
                )
                raise


if __name__ == "__main__":
    main()
