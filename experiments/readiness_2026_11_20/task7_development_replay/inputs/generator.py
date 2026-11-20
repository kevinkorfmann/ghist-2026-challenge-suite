"""Independent small-segment biological checks, not chromosome calibration.

All parameters are fixed before running and unrelated to challenge answers.
No result from this diagnostic is a calibrated submission or matched power grid.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
import tempfile
import warnings
from pathlib import Path

import msprime
import numpy as np
import stdpopsim

from ghist_single_sweep.peaks import nearest_exonic_position, read_raisd, separated_peaks
from ghist_single_sweep.recovery import localize_from_scan
from ghist_single_sweep.scan import VariantMatrix, scan_windows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--coding-bed", type=Path, required=True)
    parser.add_argument("--slim-path", type=Path, required=True)
    parser.add_argument("--raisd", type=Path, required=True)
    parser.add_argument("--sweeps", type=int, choices=(0, 1, 2), required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--scale", type=float, default=10)
    args = parser.parse_args()
    left, right = 36_000_000, 37_000_000
    intervals = []
    for line in args.coding_bed.read_text().splitlines():
        fields = line.split()
        if fields and not line.startswith("#") and fields[0].removeprefix("chr") == "7":
            a, b = max(left, int(fields[1])), min(right, int(fields[2]))
            if a < b:
                intervals.append((a, b))
    if not intervals:
        raise ValueError("No exons in the simulation segment")
    # Biological quantities are passed directly: SLiM performs its own Q scaling.
    demography = msprime.Demography()
    demography.add_population(name="population", initial_size=12000, growth_rate=math.log(10) / 300)
    demography.add_population_parameters_change(time=300, population="population", initial_size=1200, growth_rate=0)
    demography.add_population_parameters_change(time=1200, population="population", initial_size=6000, growth_rate=0)
    model = stdpopsim.DemographicModel(
        id="IndependentBG_1", description="Independent bottleneck and recent growth",
        long_description="Arbitrary diagnostic history, without challenge answers.",
        generation_time=25, model=demography, mutation_rate=1.44e-8,
    )
    species = stdpopsim.get_species("HomSap")
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        contig = species.get_contig("7", left=left, right=right, genetic_map="HapMapII_GRCh38", mutation_rate=1.44e-8)
    if any("Genetic map not found" in str(w.message) for w in captured):
        raise ValueError("Requested genetic map missing: refusing flat-map fallback")
    contig.add_dfe(np.asarray(intervals, dtype=np.int64), species.get_dfe("Gamma_K17"))
    targets = [left + 370000] if args.sweeps == 1 else [left + 270000, left + 730000]
    positions = [nearest_exonic_position(p, intervals) for p in targets[:args.sweeps]]
    events = []
    for index, position in enumerate(positions):
        locus = f"sweep_{index}"
        contig.add_single_site(id=locus, coordinate=position)
        events.extend(stdpopsim.selective_sweep(
            single_site_id=locus, population="population", mutation_generation_ago=650 + index * 100,
            selection_coeff=0.04, min_freq_at_end=0.95,
        ))
    ts = stdpopsim.get_engine("slim").simulate(
        model, contig, samples={"population": 25}, seed=args.seed,
        extended_events=events, slim_path=str(args.slim_path),
        slim_scaling_factor=args.scale, slim_burn_in=1.0,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    ts.dump(args.out.with_suffix(".trees"))
    matrix = ts.genotype_matrix()
    sites = np.asarray([site.position - left for site in ts.sites()], dtype=np.int64)
    keep = np.all((matrix == 0) | (matrix == 1), axis=1) & (matrix.sum(axis=1) > 0) & (matrix.sum(axis=1) < matrix.shape[1])
    variants = VariantMatrix("7", right - left, sites[keep], matrix[keep].astype(np.uint8), int((~keep).sum()))
    windows = scan_windows(variants, np.asarray([[0, right - left]]), window_size=100000, step_size=20000)
    localization = localize_from_scan(variants, windows)
    valid = np.flatnonzero(np.isfinite(windows["composite_score"]))
    peak = valid[np.argmax(windows["composite_score"][valid])]
    composite = int(windows["center"][peak]) + left
    estimates = {"composite": composite, "flank_ld_100kb": localization["flank_ld_100kb"] + left,
                 "largest_gap_midpoint": localization["largest_gap_midpoint"] + left}
    estimates["median_all_three"] = float(np.median(list(estimates.values())))
    with tempfile.TemporaryDirectory(dir=args.out.parent) as temporary:
        folder = Path(temporary).resolve()
        with (folder / "simulation.vcf").open("w") as handle:
            ts.write_vcf(handle, contig_id="7", position_transform=lambda pos: np.floor(pos).astype(int) + 1)
        subprocess.run([str(args.raisd.resolve()), "-n", "scientific", "-I", str(folder / "simulation.vcf"), "-w", "50", "-R", "-f"], cwd=folder, check=True, stdout=subprocess.DEVNULL)
        report_path = args.out.with_suffix(".raisd.txt")
        report_path.write_bytes((folder / "RAiSD_Report.scientific").read_bytes())
    report = read_raisd(report_path)
    peak_rows = separated_peaks(report, 100000)[:10]
    peaks = [{"position": int(report[i, 0]) - 1, "score": float(report[i, 6])} for i in peak_rows]
    payload = {
        "status": "diagnostic only; no full-chromosome calibration", "seed": args.seed,
        "generator": {"stdpopsim": stdpopsim.__version__, "demography": "bottleneck_growth",
                      "ancestral_size": 6000, "bottleneck_size": 1200, "present_size": 12000,
                      "bottleneck_start": 1200, "growth_start": 300, "dfe": "Gamma_K17",
                      "background_selection": True, "genetic_map": "HapMapII_GRCh38",
                      "mutation_rate": 1.44e-8, "slim_scaling_factor": args.scale,
                      "slim_burn_in": 1.0, "left": left, "right": right, "diploid_samples": 25,
                      "coding_bed_sha256": hashlib.sha256(args.coding_bed.read_bytes()).hexdigest()},
        "truth_positions": positions, "estimates": estimates,
        "localization": localization["localization"], "raisd_peaks": peaks,
        "single_sweep_errors": {key: abs(value - positions[0]) for key, value in estimates.items()} if len(positions) == 1 else None,
        "accepted_sites": len(variants.positions), "public_testing_answers_used": False,
        "submission_count": 0,
    }
    args.out.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
