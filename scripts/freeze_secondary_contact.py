from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from ghist_secondary_contact.contract import require_spectrum_contract


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


def require_resolution_gate(source: dict, label: str) -> None:
    if not source.get("selection_valid"):
        raise ValueError(f"{label} resolution comparison did not pass boundary checks")
    if source.get("selected_resolution") != "p40":
        raise ValueError(
            f"{label} resolution comparison selected "
            f"{source.get('selected_resolution')!r}, expected 'p40'"
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-fit", type=Path, required=True)
    parser.add_argument("--standard-fit", type=Path, required=True)
    parser.add_argument("--calibration", type=Path, required=True)
    parser.add_argument("--resolution", type=Path, required=True)
    parser.add_argument("--broad-resolution", type=Path, required=True)
    parser.add_argument("--varying-resolution", type=Path, required=True)
    parser.add_argument("--stress-comparison", type=Path, required=True)
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
        "varying_resolution": args.varying_resolution,
        "stress_comparison": args.stress_comparison,
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
    varying_resolution = sources["varying_resolution"]
    stress_comparison = sources["stress_comparison"]
    transfer = sources["stress_transfer"]
    bootstrap = sources["bootstrap"]
    audit = sources["vcf_audit"]
    for evidence in (
        calibration,
        resolution,
        broad_resolution,
        varying_resolution,
        stress_comparison,
        transfer,
        bootstrap,
    ):
        require_spectrum_contract(evidence)

    standard_values = selected_values(standard)
    refined_values = selected_values(raw)
    numerical_relative_change = {
        key: refined_values[key] / standard_values[key] - 1.0 for key in PARAMETER_KEYS
    }
    prediction = calibration["observed_accepted"]
    for label, comparison in (
        ("fixed-center", resolution),
        ("broad", broad_resolution),
        ("varying-demography", varying_resolution),
    ):
        require_resolution_gate(comparison, label)
    if not calibration["untouched_evaluation"].get("correction_accepted"):
        raise ValueError(
            "point-grouped transfer correction did not pass its acceptance gate"
        )
    if not stress_comparison.get("selection_valid"):
        raise ValueError("raw stress comparison contains optimizer-boundary results")
    if max(abs(value) for value in numerical_relative_change.values()) >= 0.01:
        raise ValueError(
            "refined integration changed at least one parameter by 1% or more"
        )
    nominal = transfer.get("conditions", {}).get("nominal")
    if nominal is None:
        raise ValueError("stress transfer is missing the nominal condition")
    if (
        nominal["corrected"]["mean_competition_rrmse"]
        >= nominal["raw"]["mean_competition_rrmse"]
    ):
        raise ValueError("correction did not improve mean nominal stress risk")
    result = {
        "challenge": "Secondary Contact",
        "status": "method frozen; no submission",
        "submission_count": 0,
        "selected_method": {
            "spectrum": "unfolded joint frequency spectrum under historical hidden-epoch convention; current generator match unverified",
            "projection": raw["projection"],
            "integration_dt_fac": raw["integration_dt_fac"],
            "adaptive_integration": raw["adaptive_integration"],
            "finite_sites": "JC69 mutations on a discrete genome with multiallelic filtering",
            "correction": calibration["accepted_prediction_method"],
            "correction_anchors": calibration["extra_calibration_patterns"],
        },
        "testing_prediction": prediction,
        "raw_refined_fit": refined_values,
        "conservative_sampling_interval_95": bootstrap["observed_sampling_interval_95"],
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
            "varying_demography_resolution": {
                "selected": varying_resolution["selected_resolution"],
                "replicates": len(varying_resolution["replicates"]),
                "summaries": varying_resolution["summaries"],
                "paired_comparisons": varying_resolution["paired_comparisons"],
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
            "raw_stress": stress_comparison["summaries"],
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
