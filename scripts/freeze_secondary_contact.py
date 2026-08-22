from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


PARAMETER_KEYS = (
    "population1_effective_size",
    "population2_effective_size",
    "generations_since_split",
    "generations_since_contact",
    "migration_rate",
)


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def selected_values(source: dict) -> dict:
    return {key: source[key] for key in PARAMETER_KEYS}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-fit", type=Path, required=True)
    parser.add_argument("--standard-fit", type=Path, required=True)
    parser.add_argument("--calibration", type=Path, required=True)
    parser.add_argument("--resolution", type=Path, required=True)
    parser.add_argument("--broad-resolution", type=Path, required=True)
    parser.add_argument("--stress-transfer", type=Path, required=True)
    parser.add_argument("--bootstrap", type=Path, required=True)
    parser.add_argument("--vcf-audit", type=Path, required=True)
    parser.add_argument("--job-id", action="append", default=[])
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    paths = {
        "raw_fit": args.raw_fit,
        "standard_fit": args.standard_fit,
        "calibration": args.calibration,
        "resolution": args.resolution,
        "broad_resolution": args.broad_resolution,
        "stress_transfer": args.stress_transfer,
        "bootstrap": args.bootstrap,
        "vcf_audit": args.vcf_audit,
    }
    sources = {name: load(path) for name, path in paths.items()}
    raw = sources["raw_fit"]
    standard = sources["standard_fit"]
    calibration = sources["calibration"]
    resolution = sources["resolution"]
    broad_resolution = sources["broad_resolution"]
    transfer = sources["stress_transfer"]
    bootstrap = sources["bootstrap"]
    audit = sources["vcf_audit"]

    standard_values = selected_values(standard)
    refined_values = selected_values(raw)
    numerical_relative_change = {
        key: refined_values[key] / standard_values[key] - 1.0
        for key in PARAMETER_KEYS
    }
    prediction = calibration["observed_accepted"]
    result = {
        "challenge": "Secondary Contact",
        "status": "method frozen; no submission",
        "submission_count": 0,
        "selected_method": {
            "spectrum": "generator-matched unfolded joint frequency spectrum",
            "projection": raw["projection"],
            "integration_dt_fac": raw["integration_dt_fac"],
            "adaptive_integration": raw["adaptive_integration"],
            "finite_sites": "JC69 mutations on a discrete genome with multiallelic filtering",
            "correction": calibration["accepted_prediction_method"],
            "correction_anchors": calibration["extra_calibration_patterns"],
        },
        "testing_prediction": prediction,
        "raw_refined_fit": refined_values,
        "conservative_sampling_interval_95": bootstrap[
            "observed_sampling_interval_95"
        ],
        "validation": {
            "resolution": {
                "selected": resolution["selected_resolution"],
                "replicates": len(resolution["replicates"]),
                "summaries": resolution["summaries"],
                "paired_comparisons": resolution["paired_comparisons"],
            },
            "broad_resolution": {
                "replicates": len(broad_resolution["replicates"]),
                "ranking": broad_resolution["ranking"],
                "summaries": broad_resolution["summaries"],
            },
            "point_grouped_transfer": {
                "partition_point_counts": calibration["point_grouped_partitions"],
                "selected_method": calibration["selected_method"],
                "untouched_evaluation": calibration["untouched_evaluation"],
            },
            "numerical_convergence": {
                "standard_fit": standard_values,
                "refined_fit": refined_values,
                "refined_minus_standard_relative": numerical_relative_change,
            },
            "stress_transfer": transfer["conditions"],
        },
        "data_provenance": {
            "testing_vcf_sha256": audit["population"]["file_sha256"],
            "retained_biallelic_records": audit["population"]["records"],
            "filtered_multiallelic_records": audit["population"][
                "skipped_non_biallelic_records"
            ],
            "ancestral_reference_records": audit["population"][
                "ancestral_reference_records"
            ],
        },
        "betty_job_ids": args.job_id,
        "source_artifacts": {
            name: {"path": str(path), "sha256": sha256(path)}
            for name, path in paths.items()
        },
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
