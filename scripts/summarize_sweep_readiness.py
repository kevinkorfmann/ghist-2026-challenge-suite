"""Summarize frozen sweep simulations, retaining failures and missing cases."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from pathlib import Path

import numpy as np
from scipy.stats import beta

from ghist_submission.validate import validate_submission
from ghist_single_sweep.peaks import highest_peak_position, read_raisd


def binomial_interval(successes: int, total: int) -> list[float] | None:
    if not total:
        return None
    return [
        0.0
        if successes == 0
        else float(beta.ppf(0.025, successes, total - successes + 1)),
        1.0
        if successes == total
        else float(beta.ppf(0.975, successes + 1, total - successes)),
    ]


def metrics(records):
    if not records:
        return {"completed": 0}
    errors = np.asarray([r["absolute_error"] for r in records])
    references = np.asarray([r["reference_absolute_error"] for r in records])
    eligible = [r for r in records if r["consensus_eligible"]]
    covered = sum(r["covered"] for r in records)
    eligible_covered = sum(r["covered"] for r in eligible)
    # This bound does not assume equal success probabilities across fixed strata.
    radius = math.sqrt(math.log(2 / 0.05) / (2 * len(records)))
    rate = covered / len(records)
    return {
        "completed": len(records),
        "covered": covered,
        "coverage": rate,
        "coverage_95pct_hoeffding_independent_bounded": [
            max(0, rate - radius),
            min(1, rate + radius),
        ],
        "coverage_95pct_clopper_pearson_if_exchangeable": binomial_interval(
            covered, len(records)
        ),
        "median_absolute_error": float(np.median(errors)),
        "q90_absolute_error": float(np.quantile(errors, 0.9)),
        "maximum_absolute_error": float(np.max(errors)),
        "consensus_eligible": len(eligible),
        "consensus_covered": eligible_covered,
        "consensus_conditional_coverage": eligible_covered / len(eligible)
        if eligible
        else None,
        "consensus_coverage_95pct_clopper_pearson_if_exchangeable": binomial_interval(
            eligible_covered, len(eligible)
        ),
        "abstention_fraction": 1 - len(eligible) / len(records),
        "reference_median_absolute_error": float(np.median(references)),
        "reference_q90_absolute_error": float(np.quantile(references, 0.9)),
        "reference_within_60kb": int(sum(references < 60_000)),
        "valid_formats": sum(r["format_validation"]["format_valid"] for r in records),
    }


def summarize(design_path: Path, directory: Path) -> dict:
    design = json.loads(design_path.read_text())
    cases = design["cases"]
    if (
        not cases
        or len({c["case"] for c in cases}) != len(cases)
        or len({c["seed"] for c in cases}) != len(cases)
    ):
        raise ValueError("Design requires distinct nonempty case identifiers and seeds")
    digest = hashlib.sha256(design_path.read_bytes()).hexdigest()
    records, statuses = [], []
    for case in cases:
        folder = directory / f"case_{case['case']:03d}"
        status_path = folder / "status.json"
        if not status_path.exists():
            statuses.append({"case": case["case"], "state": "NOT_RUN"})
            continue
        status = json.loads(status_path.read_text())
        if (
            status["design_sha256"] != digest
            or status["seed"] != case["seed"]
            or status["case"] != case["case"]
        ):
            raise ValueError("Case status does not match the frozen design")
        if status["state"] not in {"COMPLETED", "FAILED", "RUNNING"}:
            raise ValueError("Unknown case execution state")
        statuses.append(status)
        if status["state"] == "COMPLETED":
            record = json.loads((folder / "metrics.json").read_text())
            for key in (
                "case",
                "seed",
                "history",
                "selection_coefficient",
                "origin_generations",
            ):
                if record[key] != case[key]:
                    raise ValueError(f"Metrics identity mismatch: {key}")
            if record["method"] != design["predictor"]["method"]:
                raise ValueError("Prediction method differs from frozen design")
            prediction_path, bed_path = (
                folder / "prediction.json",
                folder / "candidate.bed",
            )
            if (
                hashlib.sha256(prediction_path.read_bytes()).hexdigest()
                != record["prediction_sha256"]
                or hashlib.sha256(bed_path.read_bytes()).hexdigest()
                != record["bed_sha256"]
            ):
                raise ValueError("Prediction artifact hash mismatch")
            prediction = json.loads(prediction_path.read_text())
            validation = validate_submission(bed_path, 5, "testing")
            if (
                prediction["method"] != record["method"]
                or prediction["localization"]["truth_used_for_inference"] is not False
            ):
                raise ValueError("Prediction provenance mismatch")
            lo, hi = prediction["interval"]
            if validation["intervals"] != [(lo, hi)]:
                raise ValueError("BED and prediction intervals differ")
            reference_hashes = {}
            for label, filename in (
                ("report", "RAiSD_Report.reference"),
                ("info", "RAiSD_Info.reference"),
            ):
                key = f"reference_{label}_sha256"
                reference_hashes[key] = hashlib.sha256(
                    (folder / filename).read_bytes()
                ).hexdigest()
                if key in record and reference_hashes[key] != record[key]:
                    raise ValueError("Reference artifact hash mismatch")
            reference_position = highest_peak_position(
                read_raisd(folder / "RAiSD_Report.reference"),
                prediction["contig_length"],
            )
            reference_sites = [
                int(value)
                for value in re.findall(
                    r"\|\s+snps\s+(\d+)", (folder / "RAiSD_Info.reference").read_text()
                )
            ]
            if reference_sites != [prediction["accepted_biallelic_sites"]]:
                raise ValueError("Reference and prediction SNP counts differ")
            expected = {
                "absolute_error": abs(prediction["center"] - case["position"]),
                "reference_position": reference_position,
                "reference_absolute_error": abs(reference_position - case["position"]),
                "covered": lo <= case["position"] < hi,
                "interval_width": hi - lo,
                "consensus_eligible": prediction["consensus_eligible"],
            }
            if any(record[key] != value for key, value in expected.items()):
                raise ValueError(
                    "Stored metrics differ from written predictions and frozen truth"
                )
            record["format_validation"] = validation
            record["reference_hash_attestation"] = (
                "driver and summary"
                if "reference_report_sha256" in record
                else "summary reread of archived raw report/info; legacy driver omitted hashes"
            )
            record.update(reference_hashes)
            records.append(record)
    strata = {}
    for case in design["cases"]:
        key = f"{case['history']}_s{case['selection_coefficient']}"
        if key not in strata:
            strata[key] = metrics(
                [
                    r
                    for r in records
                    if r["history"] == case["history"]
                    and r["selection_coefficient"] == case["selection_coefficient"]
                ]
            )
    return {
        "status": "simulation evidence only; no final prediction or scientific certification",
        "design_sha256": digest,
        "planned": len(design["cases"]),
        "completed": len(records),
        "operational_failures": sum(s["state"] == "FAILED" for s in statuses),
        "not_run": sum(s["state"] == "NOT_RUN" for s in statuses),
        "all_planned_complete": len(records) == len(design["cases"]),
        "completed_only_metrics": metrics(records),
        "completed_only_scope": "Accuracy is conditional on successful execution. Planned-denominator quantities below prevent silently dropping failures; unrun/running cases are unknown, not observed biological failures.",
        "planned_denominator": {
            "planned": len(cases),
            "completed": len(records),
            "failed": sum(s["state"] == "FAILED" for s in statuses),
            "running": sum(s["state"] == "RUNNING" for s in statuses),
            "not_run": sum(s["state"] == "NOT_RUN" for s in statuses),
            "completed_valid_formats": sum(
                r["format_validation"]["format_valid"] for r in records
            ),
            "valid_format_yield_lower_bound": sum(
                r["format_validation"]["format_valid"] for r in records
            )
            / len(cases),
            "covered_completed_cases": sum(r["covered"] for r in records),
            "coverage_success_yield_lower_bound": sum(r["covered"] for r in records)
            / len(cases),
            "interpretation": "Known successes divided by all planned cases; incomplete lower bounds, not estimates that unrun cases are biological misses. Operational failures remain in the planned denominator.",
        },
        "strata": strata,
        "uncertainty_scope": "Clopper-Pearson is exact for exchangeable Bernoulli replicates within a parameter stratum. Pooled heterogeneous strata additionally receive a conservative Hoeffding bound; descriptive conditional intervals do not certify the selection rule.",
        "statuses": statuses,
        "records": records,
        "submission_count": 0,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--design", type=Path, required=True)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = summarize(args.design, args.input_dir)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps(
            {
                key: value
                for key, value in result.items()
                if key not in {"statuses", "records"}
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
