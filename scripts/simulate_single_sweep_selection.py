from __future__ import annotations

import argparse
import csv
import json
import warnings
from pathlib import Path

import numpy as np
import stdpopsim

from ghist_single_sweep.scan import VariantMatrix, scan_windows
from ghist_single_sweep.recovery import localize_from_scan


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
    parser.add_argument("--selection-start-generation", type=float)
    parser.add_argument("--minimum-start-frequency", type=float)
    parser.add_argument("--selection-coefficient", type=float, default=0.05)
    parser.add_argument("--minimum-end-frequency", type=float, default=0.95)
    parser.add_argument("--samples", type=int, default=25)
    parser.add_argument(
        "--mutation-rate",
        type=float,
        default=1.44e-8,
        help="Per generation rate; testing rate by default, set 1.29e-8 for final",
    )
    parser.add_argument("--slim-scaling-factor", type=float, default=1.0)
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
            mutation_rate=args.mutation_rate,
            recombination_rate=species.genome.mean_recombination_rate,
        )
    else:
        with warnings.catch_warnings(record=True) as captured:
            warnings.simplefilter("always")
            contig = species.get_contig(
                args.chromosome,
                genetic_map=args.genetic_map,
                left=args.left,
                right=args.left + args.length,
                mutation_rate=args.mutation_rate,
            )
        if any("Genetic map not found" in str(w.message) for w in captured):
            raise ValueError(
                "Requested genetic map missing: refusing flat-map fallback"
            )
    locus_id = "single_hard_sweep"
    contig.add_single_site(id=locus_id, coordinate=args.position)
    sweep_options = {
        "single_site_id": locus_id,
        "population": "generic",
        "mutation_generation_ago": args.mutation_generation,
        "selection_coeff": args.selection_coefficient,
        "min_freq_at_end": args.minimum_end_frequency,
    }
    if args.selection_start_generation is not None:
        sweep_options["start_generation_ago"] = args.selection_start_generation
    if args.minimum_start_frequency is not None:
        sweep_options["min_freq_at_start"] = args.minimum_start_frequency
    extended_events = stdpopsim.selective_sweep(
        **sweep_options,
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
    composite_index = np.flatnonzero(valid)[
        np.argmax(windows["composite_score"][valid])
    ]
    soft_valid = np.isfinite(windows["soft_composite_score"])
    soft_index = np.flatnonzero(soft_valid)[
        np.argmax(windows["soft_composite_score"][soft_valid])
    ]
    h12_valid = np.isfinite(windows["h12"])
    h12_index = np.flatnonzero(h12_valid)[np.argmax(windows["h12"][h12_valid])]
    h12_h2_score = windows["z_high_h12"] + 0.5 * np.maximum(
        windows["z_high_h2_h1"], 0.0
    )
    joint_valid = np.isfinite(h12_h2_score)
    joint_index = np.flatnonzero(joint_valid)[np.argmax(h12_h2_score[joint_valid])]
    localized = localize_from_scan(variants, windows)
    estimates = {
        "composite": int(windows["center"][composite_index] + args.left),
        "soft_composite": int(windows["center"][soft_index] + args.left),
        "h12": int(windows["center"][h12_index] + args.left),
        "h12_h2": int(windows["center"][joint_index] + args.left),
        "flank_ld_100kb": int(localized["flank_ld_100kb"] + args.left),
        "largest_gap_midpoint": localized["largest_gap_midpoint"] + args.left,
    }
    estimate_values = np.asarray(list(estimates.values()), dtype=float)
    median_pi = float(np.nanmedian(windows["pi"]))
    payload = {
        "state": "COMPLETED",
        "localization": localized["localization"],
        "background_selection": False,
        "demography": "Zigzag_1S14",
        "mutation_rate": contig.mutation_rate,
        "stdpopsim_version": stdpopsim.__version__,
        "truth": {
            "position": args.position,
            "mutation_generation": args.mutation_generation,
            "selection_coefficient": args.selection_coefficient,
            "minimum_end_frequency": args.minimum_end_frequency,
            "selection_start_generation": args.selection_start_generation,
            "minimum_start_frequency": args.minimum_start_frequency,
        },
        "estimates": estimates,
        "absolute_error": {
            name: abs(float(value) - args.position) for name, value in estimates.items()
        },
        "diagnostics": {
            "estimator_spread": float(np.ptp(estimate_values)),
            "largest_gap_length": localized["largest_gap_length"],
            "peak_composite_score": float(windows["composite_score"][composite_index]),
            "peak_h12": float(windows["h12"][composite_index]),
            "maximum_h12": float(np.nanmax(windows["h12"])),
            "soft_peak_h12": float(windows["h12"][soft_index]),
            "soft_peak_h2_h1": float(windows["h2_h1"][soft_index]),
            "peak_soft_composite_score": float(
                windows["soft_composite_score"][soft_index]
            ),
            "peak_pi_relative_median": float(
                windows["pi"][composite_index] / median_pi
            ),
            "peak_tajima_d": float(windows["tajima_d"][composite_index]),
        },
        "accepted_biallelic_sites": int(variants.positions.size),
        "skipped_sites": variants.skipped_records,
        "sequence_length": args.length,
        "scan_window_size": 100_000,
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
