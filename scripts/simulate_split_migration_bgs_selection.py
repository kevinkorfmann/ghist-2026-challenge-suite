from __future__ import annotations

import argparse
import json
from pathlib import Path

import msprime
import numpy as np
import stdpopsim

from ghist_growth_bgs.audit import coding_distance, distance_label, load_coding_intervals
from ghist_growth_bgs.blocks import distance_stratum_segments
from ghist_split_migration_bgs.spectrum import fit_isolation_with_migration_spectrum


MIDDLE_LABELS = {"noncoding_1kb_10kb", "noncoding_10kb_100kb"}


def overlap_length(
    intervals: list[tuple[int, int]], start: int, end: int
) -> int:
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
    args = parser.parse_args()
    if not 0 <= args.left < args.right:
        raise ValueError("simulation coordinates are invalid")
    if min(
        args.ancestral_size,
        args.afr_size,
        args.eur_size,
        args.generations,
        args.migration_rate,
    ) <= 0:
        raise ValueError("demographic parameters must be positive")

    scale = args.slim_scaling_factor
    mutation_unscaled = args.mutation_rate / scale
    demography = msprime.Demography()
    demography.add_population(name="ANC", initial_size=args.ancestral_size * scale)
    demography.add_population(name="AFR", initial_size=args.afr_size * scale)
    demography.add_population(name="EUR", initial_size=args.eur_size * scale)
    demography.set_symmetric_migration_rate(
        ["AFR", "EUR"], rate=args.migration_rate / scale
    )
    demography.add_population_split(
        time=args.generations * scale,
        derived=["AFR", "EUR"],
        ancestral="ANC",
    )
    model = stdpopsim.DemographicModel(
        id="SplitMigration_2P",
        description="Two populations with continuous symmetric migration",
        long_description="A constant-size ancestor splits into two constant-size populations.",
        generation_time=0.1,
        model=demography,
        mutation_rate=mutation_unscaled,
    )
    species = stdpopsim.get_species("DroMel")
    contig = species.get_contig(
        args.chromosome,
        genetic_map="ComeronCrossoverV2_dm6",
        left=args.left,
        right=args.right,
        mutation_rate=mutation_unscaled,
    )
    scale_recombination_map(contig, 0.5)
    coding = load_coding_intervals(args.coding_bed)
    contig.add_dfe(
        np.asarray(coding[args.chromosome], dtype=np.int64),
        species.get_dfe("Gamma_H17"),
    )
    engine = stdpopsim.get_engine("slim")
    tree_sequence = engine.simulate(
        model,
        contig,
        samples={"AFR": args.afr_samples, "EUR": args.eur_samples},
        seed=args.seed,
        slim_path=str(args.slim_path),
        slim_scaling_factor=scale,
        slim_burn_in=args.slim_burn_in,
    )
    args.trees.parent.mkdir(parents=True, exist_ok=True)
    tree_sequence.dump(args.trees)

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
    fit = fit_isolation_with_migration_spectrum(
        spectrum,
        mutation_rate=args.mutation_rate,
        sequence_length=exposure,
        projection=(60, 50),
        restarts=args.fit_restarts,
        workers=1,
        maxiter=600,
        seed=args.seed + 1,
    )
    truth = {
        "ancestral_population_size": args.ancestral_size,
        "AFR_population_size": args.afr_size,
        "EUR_population_size": args.eur_size,
        "generations_since_split": args.generations,
        "migration_rate": args.migration_rate,
    }
    payload = {
        "truth_in_analysis_units": truth,
        "generator": {
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
        "relative_error": {
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
    args.out.write_text(json.dumps(payload, indent=2) + "\n")
    np.savez_compressed(args.out.with_suffix(".sfs.npz"), spectrum=spectrum)
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
