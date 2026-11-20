from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

from ghist_inference.bootstrap import require_summary
from ghist_inference.checkpoint import atomic_json
from ghist_inference.freeze import (
    continuous_prediction,
    require_optimizer_convergence,
    verified_projection_median,
    validate_physical_identity,
)
from ghist_growth_bgs.simulation import require_physical_scaling


PARAMETERS = (
    "AFR_population_size",
    "EUR_population_size",
    "generations_since_split",
    "migration_rate",
)
PROJECTIONS = ((50, 40), (70, 50), (94, 68))


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--refined", type=Path, required=True)
    parser.add_argument("--bootstrap", type=Path, required=True)
    parser.add_argument("--calibration", type=Path, required=True)
    parser.add_argument("--near-recovery", type=Path, required=True)
    parser.add_argument("--distance", type=Path, action="append", default=[])
    parser.add_argument("--arm", type=Path, action="append", default=[])
    parser.add_argument("--figure", type=Path, action="append", default=[])
    parser.add_argument("--job-id", action="append", default=[])
    parser.add_argument("--testing-vcf-sha256", required=True)
    parser.add_argument("--coding-bed-sha256", required=True)
    parser.add_argument("--recombination-map-sha256", required=True)
    parser.add_argument("--final-sealed-sha256", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    baseline = load(args.baseline)
    refined = load(args.refined)
    bootstrap = load(args.bootstrap)
    calibration = load(args.calibration)
    recovery = load(args.near_recovery)
    require_physical_scaling([calibration, recovery])
    raw = verified_projection_median(
        baseline.get("by_projection"), baseline.get("projection_median"), PARAMETERS,
        label="split fit projections",
        expected_projections=[f"{left}x{right}" for left, right in PROJECTIONS],
    )
    for fit in baseline["by_projection"].values():
        validate_physical_identity(fit, 4)
    require_optimizer_convergence(refined.get("fit", refined), "refined split fit")
    validate_physical_identity(refined.get("fit", refined), 4)

    if calibration.get("selected_model") != "identity":
        raise ValueError("held-out calibration selected a correction")
    limits = {
        "AFR_population_size": 0.05,
        "EUR_population_size": 0.05,
        "generations_since_split": 0.05,
        "migration_rate": 0.15,
    }
    for parameter, limit in limits.items():
        if (
            not math.isfinite(recovery["relative_error"][parameter])
            or abs(recovery["relative_error"][parameter]) >= limit
        ):
            raise ValueError(f"{parameter} matched recovery exceeds tolerance")

    refined_fit = continuous_prediction(refined.get("fit", refined), PARAMETERS)
    reference = continuous_prediction(baseline["by_projection"]["50x40"], PARAMETERS)
    numerical_change = {
        parameter: refined_fit[parameter] / reference[parameter] - 1.0
        for parameter in PARAMETERS
    }
    if max(abs(value) for value in numerical_change.values()) >= 0.01:
        raise ValueError("refined integration changes a parameter by at least 1%")

    interval = require_summary(bootstrap, 4, observed_fits=baseline["by_projection"])
    source_paths = {
        "baseline": args.baseline,
        "refined": args.refined,
        "bootstrap": args.bootstrap,
        "calibration": args.calibration,
        "near_recovery": args.near_recovery,
        **{f"distance_{i}": path for i, path in enumerate(args.distance)},
        **{f"arm_{i}": path for i, path in enumerate(args.arm)},
    }
    result = {
        "challenge": "Split migration with background selection",
        "status": "method frozen; no submission",
        "submission_count": 0,
        "selected_method": {
            "model": "two-population isolation with continuous symmetric migration",
            "spectrum": "unfolded 1–100 kb noncoding joint frequency spectrum",
            "projections": [list(projection) for projection in PROJECTIONS],
            "aggregation": "continuous projection median; no preset rounding",
            "calibration": "identity; no correction",
        },
        "testing_prediction": raw,
        "raw_projection_median": {key: raw[key] for key in PARAMETERS},
        "linked_block_interval_95": interval,
        "validation": {
            "projection_fits": baseline["by_projection"],
            "matched_selection_recovery": recovery,
            "held_out_calibration": {
                "validation": calibration["validation"],
                "truth_points": calibration["truth_points"],
                "selected_model": calibration["selected_model"],
                "evaluations": calibration["evaluations"],
            },
            "linked_block_replicates": bootstrap["replicates"],
            "bootstrap_contract": bootstrap["bootstrap_contract"],
            "bootstrap_planned_replicates": bootstrap["planned_replicates"],
            "interval_coverage_validated": False,
            "chromosome_arm_fits": [load(path) for path in args.arm],
            "numerical_refinement_relative_change": numerical_change,
        },
        "data_provenance": {
            "testing_vcf_sha256": args.testing_vcf_sha256,
            "coding_bed_sha256": args.coding_bed_sha256,
            "recombination_map_sha256": args.recombination_map_sha256,
            "final_vcf_sealed_sha256": args.final_sealed_sha256,
            "segregating_biallelic_sites": baseline["metadata"][
                "ancestral_reference_records"
            ],
            "haplotypes": baseline["metadata"]["haplotypes"],
            "middle_distance_exposure_bp": baseline["physical_exposure_bp"],
        },
        "betty_job_ids": args.job_id,
        "source_artifacts": {
            name: {"path": str(path), "sha256": sha256(path)}
            for name, path in source_paths.items()
        },
        "figures": [
            {"path": str(path), "sha256": sha256(path)} for path in args.figure
        ],
    }
    atomic_json(args.out, result)
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
