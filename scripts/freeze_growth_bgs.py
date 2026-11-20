from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

from ghist_inference.freeze import continuous_prediction


PARAMETERS = ("generations", "post_growth_scaling")


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
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--middle-fit", type=Path, required=True)
    parser.add_argument("--bootstrap", type=Path, required=True)
    parser.add_argument("--convention", type=Path, required=True)
    parser.add_argument(
        "--model-comparison",
        type=Path,
        required=True,
        help="Chromosome holdout comparison from compare_growth_models.py",
    )
    parser.add_argument("--smoke", type=Path, required=True)
    parser.add_argument("--grid", type=Path, required=True)
    parser.add_argument("--calibration", type=Path, required=True)
    parser.add_argument("--stress", type=Path, required=True)
    parser.add_argument("--figure", type=Path, action="append", default=[])
    parser.add_argument("--job-id", action="append", default=[])
    parser.add_argument("--final-sealed-sha256", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    paths = {
        "audit": args.audit,
        "middle_fit": args.middle_fit,
        "bootstrap": args.bootstrap,
        "convention": args.convention,
        "model_comparison": args.model_comparison,
        "smoke": args.smoke,
        "grid": args.grid,
        "calibration": args.calibration,
        "stress": args.stress,
    }
    sources = {name: load(path) for name, path in paths.items()}
    middle = sources["middle_fit"]["strata"]["noncoding_1kb_100kb"]
    bootstrap = sources["bootstrap"]
    smoke = sources["smoke"]
    grid = sources["grid"]
    calibration = sources["calibration"]
    stress = sources["stress"]

    observed = continuous_prediction(middle["median"], PARAMETERS)
    comparison = sources["model_comparison"]
    # All legacy middle-distance, bootstrap, and SLiM calibration artifacts fit
    # exponential growth. They cannot certify a different selected convention.
    if comparison.get("selected_convention") != "exponential" or not comparison.get(
        "selection_consistent_across_groups"
    ):
        raise ValueError(
            "Chromosome holdouts do not consistently support the exponential model; regenerate fits and BGS recovery for the selected model"
        )
    if smoke.get("preferred") != "exponential":
        raise ValueError("Independent recovery does not support the fitted convention")
    if calibration.get("selected_model") != "identity":
        raise ValueError(
            "Held-out calibration selected a correction; apply and validate it before freezing"
        )
    # Check the complete recovery grid, not one handpicked central point.
    if not grid.get("results"):
        raise ValueError("Independent recovery grid is empty")
    for point, recovery in grid["results"].items():
        for parameter, limit in (("generations", 0.10), ("post_growth_scaling", 0.05)):
            error = float(recovery["relative_error"][parameter])
            if not math.isfinite(error) or abs(error) >= limit:
                raise ValueError(
                    f"Independent grid recovery exceeds tolerance: {point}, {parameter}"
                )
    if not stress.get("freeze_gate_passed"):
        raise ValueError("selection-strength and burn-in stress gate did not pass")

    interval = {
        "generations": [
            bootstrap["generations"]["quantiles"]["0.025"],
            bootstrap["generations"]["quantiles"]["0.975"],
        ],
        "post_growth_scaling": [
            bootstrap["post_growth_scaling"]["quantiles"]["0.025"],
            bootstrap["post_growth_scaling"]["quantiles"]["0.975"],
        ],
    }
    audit = sources["audit"]
    result = {
        "challenge": "Growth with background selection",
        "status": "method frozen; no submission",
        "submission_count": 0,
        "selected_method": {
            "model": "exponential growth",
            "spectrum": "unfolded middle-distance noncoding frequency spectrum",
            "distance_from_coding_bp": [1000, 100000],
            "projections": [84, 60, 40],
            "aggregation": "continuous projection median; no preset rounding",
            "calibration": "identity; no correction",
        },
        "testing_prediction": observed,
        "raw_projection_median": {
            "generations": observed["generations"],
            "post_growth_scaling": observed["post_growth_scaling"],
        },
        "linked_block_interval_95": interval,
        "validation": {
            "projection_fits": middle["by_projection"],
            "linked_block_replicates": bootstrap["replicates"],
            "known_truth_smoke": smoke,
            "all_grid_recovery": grid["results"],
            "chromosome_holdout_comparison": comparison,
            "held_out_calibration": {
                "validation": calibration["validation"],
                "selected_model": calibration["selected_model"],
                "evaluations": calibration["evaluations"],
            },
            "selection_strength_and_burnin_stress": stress,
        },
        "data_provenance": {
            "testing_vcf_sha256": audit["files"]["vcf"]["sha256"],
            "coding_bed_sha256": audit["files"]["coding_bed"]["sha256"],
            "recombination_map_sha256": audit["files"]["recombination_map"]["sha256"],
            "final_vcf_sealed_sha256": args.final_sealed_sha256,
            "segregating_biallelic_sites": audit["segregating_biallelic_sites"],
            "haplotypes": audit["expected_haplotypes"],
            "middle_distance_exposure_bp": sum(
                audit["physical_bp_by_distance_stratum"][key]
                for key in ("noncoding_1kb_10kb", "noncoding_10kb_100kb")
            ),
        },
        "betty_job_ids": args.job_id,
        "source_artifacts": {
            name: {"path": str(path), "sha256": sha256(path)}
            for name, path in paths.items()
        },
        "figures": [
            {"path": str(path), "sha256": sha256(path)} for path in args.figure
        ],
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
