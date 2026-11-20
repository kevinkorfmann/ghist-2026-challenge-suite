from __future__ import annotations

import argparse
import json
from pathlib import Path

import msprime
import numpy as np
import stdpopsim

from ghist_inference.checkpoint import (
    atomic_json,
    fit_checkpoint,
    save_spectrum_checkpoint,
    simulation_identity,
    simulation_tree_checkpoint,
)

from ghist_growth_bgs.audit import (
    coding_distance,
    distance_label,
    load_coding_intervals,
)
from ghist_growth_bgs.blocks import distance_stratum_segments
from ghist_split_migration_bgs.spectrum import fit_isolation_with_migration_spectrum
from ghist_split_migration_bgs.simulation import split_migration_demography
from ghist_growth_bgs.simulation import SIMULATION_SCALING_SCHEME


MIDDLE_LABELS = {"noncoding_1kb_10kb", "noncoding_10kb_100kb"}


def overlap_length(intervals: list[tuple[int, int]], start: int, end: int) -> int:
    return sum(
        max(0, min(right, end) - max(left, start))
        for left, right in intervals
        if left < end and right > start
    )


def scale_recombination_map(contig, factor: float) -> None:
    rates = np.asarray(contig.recombination_map.rate, dtype=np.float64)
    contig.recombination_map = msprime.RateMap(
        position=np.asarray(contig.recombination_map.position, dtype=np.float64),
        rate=np.where(np.isnan(rates), rates, rates * factor),
    )


def population_nodes(tree_sequence, name: str) -> np.ndarray:
    for population_id, population in enumerate(tree_sequence.populations()):
        metadata = population.metadata
        if isinstance(metadata, dict) and metadata.get("name") == name:
            return tree_sequence.samples(population=population_id)
    raise ValueError(f"population {name!r} is absent from the tree sequence")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--coding-bed", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--trees", type=Path, required=True)
    parser.add_argument("--slim-path", type=Path, required=True)
    parser.add_argument("--chromosome", default="2L")
    parser.add_argument("--left", type=int, default=5_000_000)
    parser.add_argument("--right", type=int, default=6_000_000)
    parser.add_argument("--ancestral-size", type=float, required=True)
    parser.add_argument("--afr-size", type=float, required=True)
    parser.add_argument("--eur-size", type=float, required=True)
    parser.add_argument("--generations", type=float, required=True)
    parser.add_argument("--migration-rate", type=float, required=True)
    parser.add_argument("--slim-scaling-factor", type=float, default=10.0)
    parser.add_argument("--slim-burn-in", type=float, default=0.1)
    parser.add_argument("--mutation-rate", type=float, default=1.015e-8)
    parser.add_argument("--afr-samples", type=int, default=47)
    parser.add_argument("--eur-samples", type=int, default=34)
    parser.add_argument("--seed", type=int, default=75194000)
    parser.add_argument("--fit-restarts", type=int, default=16)
    parser.add_argument("--skip-fit", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument(
        "--existing-trees-sha256",
        help="Required to import a tree predating simulation checkpoints",
    )
    args = parser.parse_args()
    if not 0 <= args.left < args.right:
        raise ValueError("simulation coordinates are invalid")
    if (
        min(
            args.ancestral_size,
            args.afr_size,
            args.eur_size,
            args.generations,
            args.migration_rate,
        )
        <= 0
    ):
        raise ValueError("demographic parameters must be positive")

    scale = args.slim_scaling_factor
    demography = split_migration_demography(
        args.ancestral_size,
        args.afr_size,
        args.eur_size,
        args.generations,
        args.migration_rate,
    )
    model = stdpopsim.DemographicModel(
        id="SplitMigration_2P",
        description="Two populations with continuous symmetric migration",
        long_description="A constant-size ancestor splits into two constant-size populations.",
        generation_time=0.1,
        model=demography,
        mutation_rate=args.mutation_rate,
    )
    species = stdpopsim.get_species("DroMel")
    contig = species.get_contig(
        args.chromosome,
        genetic_map="ComeronCrossoverV2_dm6",
        left=args.left,
        right=args.right,
        mutation_rate=args.mutation_rate,
    )
    scale_recombination_map(contig, 0.5)
    coding = load_coding_intervals(args.coding_bed)
    contig.add_dfe(
        np.asarray(coding[args.chromosome], dtype=np.int64),
        species.get_dfe("Gamma_H17"),
    )
    engine = stdpopsim.get_engine("slim")
    identity = simulation_identity(
        args,
        stdpopsim_version=stdpopsim.__version__,
        scaling_scheme=SIMULATION_SCALING_SCHEME,
    )

    def simulate():
        return engine.simulate(
            model,
            contig,
            samples={"AFR": args.afr_samples, "EUR": args.eur_samples},
            seed=args.seed,
            slim_path=str(args.slim_path),
            slim_scaling_factor=scale,
            slim_burn_in=args.slim_burn_in,
        )

    tree_sequence = simulation_tree_checkpoint(args, identity, simulate)

    afr_nodes = population_nodes(tree_sequence, "AFR")
    eur_nodes = population_nodes(tree_sequence, "EUR")
    sample_nodes = np.concatenate([afr_nodes, eur_nodes])
    spectrum = np.zeros((len(afr_nodes) + 1, len(eur_nodes) + 1), dtype=np.int64)
    starts = [interval[0] for interval in coding[args.chromosome]]
    accepted = multiallelic = 0
    for variant in tree_sequence.variants(samples=sample_nodes):
        if len(variant.alleles) != 2 or np.any(variant.genotypes < 0):
            multiallelic += 1
            continue
        position = int(variant.site.position)
        if not args.left <= position < args.right:
            continue
        label = distance_label(
            coding_distance(coding[args.chromosome], position, starts)
        )
        if label not in MIDDLE_LABELS:
            continue
        afr_derived = int(np.sum(variant.genotypes[: len(afr_nodes)]))
        eur_derived = int(np.sum(variant.genotypes[len(afr_nodes) :]))
        if 0 < afr_derived + eur_derived < len(sample_nodes):
            spectrum[afr_derived, eur_derived] += 1
            accepted += 1
    selected_segments = distance_stratum_segments(
        coding[args.chromosome], int(contig.length), MIDDLE_LABELS
    )
    exposure = overlap_length(selected_segments, args.left, args.right)
    save_spectrum_checkpoint(
        args.out.with_suffix(".sfs.npz"),
        spectrum,
        {"simulation": identity, "physical_exposure_bp": exposure},
    )

    truth = {
        "ancestral_population_size": args.ancestral_size,
        "AFR_population_size": args.afr_size,
        "EUR_population_size": args.eur_size,
        "generations_since_split": args.generations,
        "migration_rate": args.migration_rate,
        "mutation_rate": args.mutation_rate,
    }
    fit = None
    if not args.skip_fit:
        fit = fit_checkpoint(
            args.out.with_suffix(".fit.json"),
            spectrum=spectrum,
            fit=fit_isolation_with_migration_spectrum,
            resume=args.resume,
            settings={
                "mutation_rate": args.mutation_rate,
                "sequence_length": exposure,
                "projection": (min(60, len(afr_nodes)), min(50, len(eur_nodes))),
                "restarts": args.fit_restarts,
                "workers": 1,
                "maxiter": 600,
                "seed": args.seed + 1,
            },
        )

    payload = {
        "truth_in_analysis_units": truth,
        "generator": {
            "simulation_scaling_scheme": SIMULATION_SCALING_SCHEME,
            "stdpopsim": stdpopsim.__version__,
            "slim_scaling_factor": scale,
            "slim_burn_in": args.slim_burn_in,
            "chromosome": args.chromosome,
            "left": args.left,
            "right": args.right,
            "seed": args.seed,
            "dfe": "Gamma_H17",
            "genetic_map": "ComeronCrossoverV2_dm6 scaled by 0.5",
        },
        "spectrum": {
            "AFR_haplotypes": len(afr_nodes),
            "EUR_haplotypes": len(eur_nodes),
            "accepted_sites": accepted,
            "multiallelic_sites": multiallelic,
            "physical_exposure_bp": exposure,
        },
        "fit": fit,
        "relative_error": None
        if fit is None
        else {
            key: fit[key] / truth[key] - 1
            for key in (
                "AFR_population_size",
                "EUR_population_size",
                "generations_since_split",
                "migration_rate",
            )
        },
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    atomic_json(args.out, payload)
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
