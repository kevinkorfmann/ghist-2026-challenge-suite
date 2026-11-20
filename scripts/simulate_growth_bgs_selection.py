from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import msprime
import numpy as np
import stdpopsim

from ghist_growth_bgs.audit import (
    coding_distance,
    distance_label,
    load_coding_intervals,
)
from ghist_growth_bgs.blocks import distance_stratum_segments
from ghist_growth_bgs.spectrum import fit_growth_spectrum


MIDDLE_LABELS = {"noncoding_1kb_10kb", "noncoding_10kb_100kb"}


def overlap_length(intervals: list[tuple[int, int]], start: int, end: int) -> int:
    return sum(
        max(0, min(right, end) - max(left, start))
        for left, right in intervals
        if left < end and right > start
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--coding-bed", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--trees", type=Path, required=True)
    parser.add_argument("--slim-path", type=Path, required=True)
    parser.add_argument("--chromosome", default="2L")
    parser.add_argument("--left", type=int, default=5_000_000)
    parser.add_argument("--right", type=int, default=10_000_000)
    parser.add_argument("--ancestral-size", type=float, default=10_000)
    parser.add_argument("--generations", type=float, default=5_000)
    parser.add_argument("--scaling", type=float, default=2.0)
    parser.add_argument(
        "--convention", choices=("exponential", "sudden"), default="exponential"
    )
    parser.add_argument("--slim-scaling-factor", type=float, default=10.0)
    parser.add_argument("--slim-burn-in", type=float, default=0.1)
    parser.add_argument("--mutation-rate", type=float, default=3.36e-8)
    parser.add_argument("--samples", type=int, default=42)
    parser.add_argument("--seed", type=int, default=7593401)
    parser.add_argument("--fit-restarts", type=int, default=32)
    args = parser.parse_args()
    if not 0 <= args.left < args.right:
        raise ValueError("simulation coordinates are invalid")
    if args.scaling <= 1 or args.slim_scaling_factor < 1:
        raise ValueError("growth and simulation scaling factors are invalid")

    scale = args.slim_scaling_factor
    ancestral_unscaled = args.ancestral_size * scale
    modern_unscaled = ancestral_unscaled * args.scaling
    generations_unscaled = args.generations * scale
    mutation_unscaled = args.mutation_rate / scale
    growth_rate = (
        math.log(args.scaling) / generations_unscaled
        if args.convention == "exponential"
        else 0.0
    )
    demography = msprime.Demography()
    demography.add_population(
        name="population",
        initial_size=modern_unscaled,
        growth_rate=growth_rate,
    )
    demography.add_population_parameters_change(
        time=generations_unscaled,
        initial_size=ancestral_unscaled,
        growth_rate=0,
        population="population",
    )
    model = stdpopsim.DemographicModel(
        id="Growth_1G",
        description=f"Single-population {args.convention} growth",
        long_description=f"{args.convention.capitalize()} growth with a fixed ancestral epoch.",
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
    coding = load_coding_intervals(args.coding_bed)
    contig.add_dfe(
        np.asarray(coding[args.chromosome], dtype=np.int64),
        species.get_dfe("Gamma_H17"),
    )
    engine = stdpopsim.get_engine("slim")
    tree_sequence = engine.simulate(
        model,
        contig,
        samples={"population": args.samples},
        seed=args.seed,
        slim_path=str(args.slim_path),
        slim_scaling_factor=scale,
        slim_burn_in=args.slim_burn_in,
    )
    args.trees.parent.mkdir(parents=True, exist_ok=True)
    tree_sequence.dump(args.trees)

    starts = [interval[0] for interval in coding[args.chromosome]]
    spectrum = np.zeros(2 * args.samples + 1, dtype=np.int64)
    accepted = multiallelic = 0
    sample_nodes = tree_sequence.samples()
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
        derived = int(np.sum(variant.genotypes))
        if 0 < derived < len(sample_nodes):
            spectrum[derived] += 1
            accepted += 1
    selected_segments = distance_stratum_segments(
        coding[args.chromosome], int(contig.length), MIDDLE_LABELS
    )
    exposure = overlap_length(selected_segments, args.left, args.right)
    fit = fit_growth_spectrum(
        spectrum,
        sequence_length=exposure,
        mutation_rate=args.mutation_rate,
        projection=len(sample_nodes),
        restarts=args.fit_restarts,
        workers=1,
        maxiter=500,
        seed=args.seed + 1,
    )
    sudden_fit = fit_growth_spectrum(
        spectrum,
        sequence_length=exposure,
        mutation_rate=args.mutation_rate,
        projection=len(sample_nodes),
        restarts=args.fit_restarts,
        workers=1,
        maxiter=500,
        seed=args.seed + 2,
        convention="sudden",
    )
    matched_fit = fit if args.convention == "exponential" else sudden_fit
    payload = {
        "truth_in_analysis_units": {
            "ancestral_effective_size": args.ancestral_size,
            "generations": args.generations,
            "post_growth_scaling": args.scaling,
            "mutation_rate": args.mutation_rate,
        },
        "generator": {
            "convention": args.convention,
            "stdpopsim": stdpopsim.__version__,
            "slim_scaling_factor": scale,
            "slim_burn_in": args.slim_burn_in,
            "chromosome": args.chromosome,
            "left": args.left,
            "right": args.right,
            "seed": args.seed,
            "dfe": "Gamma_H17",
            "genetic_map": "ComeronCrossoverV2_dm6",
        },
        "spectrum": {
            "haplotypes": len(sample_nodes),
            "accepted_sites": accepted,
            "multiallelic_sites": multiallelic,
            "physical_exposure_bp": exposure,
        },
        "fit": fit,
        "matched_fit": matched_fit,
        "sudden_fit": sudden_fit,
        "log_likelihood_exponential_minus_sudden": (
            fit["multinomial_log_likelihood"] - sudden_fit["multinomial_log_likelihood"]
        ),
        "relative_error": {
            "generations": matched_fit["generations"] / args.generations - 1,
            "post_growth_scaling": matched_fit["post_growth_scaling"] / args.scaling
            - 1,
        },
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n")
    np.savez_compressed(args.out.with_suffix(".sfs.npz"), spectrum=spectrum)
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
