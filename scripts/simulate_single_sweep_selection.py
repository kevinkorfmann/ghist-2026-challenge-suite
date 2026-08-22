from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
import stdpopsim

from ghist_single_sweep.scan import VariantMatrix, flank_ld_scan, scan_windows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--trees", type=Path, required=True)
    parser.add_argument("--scan", type=Path, required=True)
    parser.add_argument("--slim-path", type=Path, required=True)
    parser.add_argument("--length", type=int, default=5_000_000)
    parser.add_argument("--position", type=int, default=2_500_000)
    parser.add_argument("--chromosome")
    parser.add_argument("--left", type=int, default=0)
    parser.add_argument("--genetic-map")
    parser.add_argument("--mutation-generation", type=float, default=1000.0)
    parser.add_argument("--selection-coefficient", type=float, default=0.05)
    parser.add_argument("--minimum-end-frequency", type=float, default=0.95)
    parser.add_argument("--samples", type=int, default=25)
    parser.add_argument("--slim-scaling-factor", type=float, default=10.0)
    parser.add_argument("--slim-burn-in", type=float, default=0.1)
    parser.add_argument("--seed", type=int, default=75221000)
    args = parser.parse_args()

    species = stdpopsim.get_species("HomSap")
    model = species.get_demographic_model("Zigzag_1S14")
    if args.chromosome is None:
        if args.left != 0 or args.genetic_map is not None:
            raise ValueError("left and genetic map require a named chromosome")
        contig = species.get_contig(
            length=args.length,
            mutation_rate=species.genome.mean_mutation_rate,
            recombination_rate=species.genome.mean_recombination_rate,
        )
    else:
        contig = species.get_contig(
            args.chromosome,
            genetic_map=args.genetic_map,
            left=args.left,
            right=args.left + args.length,
            mutation_rate=species.genome.mean_mutation_rate,
        )
    locus_id = "single_hard_sweep"
    contig.add_single_site(id=locus_id, coordinate=args.position)
    extended_events = stdpopsim.selective_sweep(
        single_site_id=locus_id,
        population="generic",
        mutation_generation_ago=args.mutation_generation,
        selection_coeff=args.selection_coefficient,
        min_freq_at_end=args.minimum_end_frequency,
    )
    engine = stdpopsim.get_engine("slim")
    tree_sequence = engine.simulate(
        model,
        contig,
        samples={"generic": args.samples},
        seed=args.seed,
        extended_events=extended_events,
        slim_path=str(args.slim_path),
        slim_scaling_factor=args.slim_scaling_factor,
        slim_burn_in=args.slim_burn_in,
    )
    args.trees.parent.mkdir(parents=True, exist_ok=True)
    tree_sequence.dump(args.trees)

    matrix = tree_sequence.genotype_matrix()
    positions = np.asarray(
        [site.position - args.left for site in tree_sequence.sites()], dtype=np.int32
    )
    keep = (
        np.all((matrix == 0) | (matrix == 1), axis=1)
        & (matrix.sum(axis=1) > 0)
        & (matrix.sum(axis=1) < matrix.shape[1])
    )
    variants = VariantMatrix(
        contig="simulated",
        contig_length=args.length,
        positions=positions[keep],
        haplotypes=matrix[keep].astype(np.uint8),
        skipped_records=int((~keep).sum()),
    )
    windows = scan_windows(
        variants,
        np.asarray([[0, args.length]], dtype=np.int64),
        window_size=100_000,
        step_size=20_000,
    )
    args.scan.parent.mkdir(parents=True, exist_ok=True)
    names = list(windows)
    with args.scan.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(names)
        writer.writerows(zip(*(windows[name] for name in names), strict=True))

    valid = np.isfinite(windows["composite_score"])
    composite_index = np.flatnonzero(valid)[np.argmax(windows["composite_score"][valid])]
    truth_local = args.position - args.left
    local_region = [max(0, truth_local - 750_000), min(args.length, truth_local + 750_000)]
    candidates = np.arange(local_region[0], local_region[1] + 1, 1000, dtype=np.int64)
    ld = flank_ld_scan(variants, candidates, flank_width=100_000)
    ld_valid = np.isfinite(ld["flank_ld_ratio"])
    ld_index = np.flatnonzero(ld_valid)[np.argmax(ld["flank_ld_ratio"][ld_valid])]
    local_lo = int(np.searchsorted(variants.positions, local_region[0], side="left"))
    local_hi = int(np.searchsorted(variants.positions, local_region[1], side="right"))
    local_positions = variants.positions[local_lo:local_hi]
    gaps = np.diff(local_positions.astype(np.int64))
    gap_index = int(np.argmax(gaps))
    estimates = {
        "composite": int(windows["center"][composite_index] + args.left),
        "flank_ld_100kb": int(ld["candidate_position"][ld_index] + args.left),
        "largest_gap_midpoint": float(
            (int(local_positions[gap_index]) + int(local_positions[gap_index + 1])) / 2
            + args.left
        ),
    }
    estimate_values = np.asarray(list(estimates.values()), dtype=float)
    median_pi = float(np.nanmedian(windows["pi"]))
    payload = {
        "state": "COMPLETED",
        "truth": {
            "position": args.position,
            "mutation_generation": args.mutation_generation,
            "selection_coefficient": args.selection_coefficient,
            "minimum_end_frequency": args.minimum_end_frequency,
        },
        "estimates": estimates,
        "absolute_error": {
            name: abs(float(value) - args.position) for name, value in estimates.items()
        },
        "diagnostics": {
            "estimator_spread": float(np.ptp(estimate_values)),
            "largest_gap_length": int(gaps[gap_index]),
            "peak_composite_score": float(windows["composite_score"][composite_index]),
            "peak_h12": float(windows["h12"][composite_index]),
            "maximum_h12": float(np.nanmax(windows["h12"])),
            "peak_pi_relative_median": float(windows["pi"][composite_index] / median_pi),
            "peak_tajima_d": float(windows["tajima_d"][composite_index]),
        },
        "accepted_biallelic_sites": int(variants.positions.size),
        "skipped_sites": variants.skipped_records,
        "sequence_length": args.length,
        "chromosome": args.chromosome,
        "left": args.left,
        "genetic_map": args.genetic_map,
        "seed": args.seed,
        "slim_scaling_factor": args.slim_scaling_factor,
        "slim_burn_in": args.slim_burn_in,
        "submission_count": 0,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
