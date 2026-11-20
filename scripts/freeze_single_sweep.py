from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from ghist_single_sweep.recovery import validate_blind_records


MAX_FROZEN_HALF_WIDTH = 75_000


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate(
    generic: dict,
    chromosome: dict,
    exact: dict,
    candidate: dict,
    raisd: dict,
) -> dict:
    for recovery in (generic, chromosome, exact):
        validate_blind_records(recovery.get("records", []))
    for label, recovery in [("generic", generic), ("chromosome", chromosome)]:
        if recovery["held_out_method_gate"]["selected"] != "median_all_three":
            raise ValueError(f"{label} recovery did not select the locked ensemble")
    generic_eval = generic["held_out_method_gate"]["evaluation_metrics"]["median_all_three"]
    chromosome_eval = chromosome["held_out_method_gate"]["evaluation_metrics"]["median_all_three"]
    if generic_eval["quantiles"]["0.9"] > 75_000:
        raise ValueError("generic held-out 90th-percentile error exceeds 75 kb")
    if chromosome_eval["quantiles"]["0.9"] > 75_000:
        raise ValueError("chromosome-map held-out 90th-percentile error exceeds 75 kb")
    if not exact["robustness_gate"]["passed"]:
        raise ValueError("untouched exact-scale recovery did not pass")

    selected = candidate["interval_selection"]
    start, end = selected["interval"]
    if selected["final_half_width"] > MAX_FROZEN_HALF_WIDTH:
        raise ValueError("candidate half-width exceeds the precision gate")
    outside = {
        name: value
        for name, value in candidate["observed_estimates"].items()
        if not start <= value < end
    }
    if outside:
        raise ValueError(f"observed estimators outside candidate interval: {outside}")
    peaks = raisd["separated_peaks_100kb"]
    reference_ratio = float(peaks[0]["mu"] / peaks[1]["mu"])
    if reference_ratio < 3.0:
        raise ValueError("reference scan lacks a three-fold separated-peak margin")
    if not start <= raisd["maximum"]["position"] < end:
        raise ValueError("reference maximum lies outside the candidate interval")
    maximum_hit_score = 1.0 / (1.0 + (end - start) * 1e-5 / 2.0)
    return {
        "generic_held_out": generic_eval,
        "chromosome_map_held_out": chromosome_eval,
        "exact_scale_recovery": {
            "replicates": exact["replicates"],
            "calibration_consensus_replicates": exact["calibration_consensus_replicates"],
            "pilot_consensus_replicates": exact["pilot_consensus_replicates"],
            "evaluation_consensus_replicates": exact["evaluation_consensus_replicates"],
            "selected_half_width": exact["interval_width_selection"]["selected_half_width"],
            "pilot_evaluation": exact["pilot_evaluation"],
            "untouched_evaluation": exact["untouched_evaluation"],
            "robustness_gate": exact["robustness_gate"],
        },
        "reference_first_to_second_peak_ratio": reference_ratio,
        "candidate_hit_score": maximum_hit_score,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--generic-recovery", type=Path, required=True)
    parser.add_argument("--chromosome-recovery", type=Path, required=True)
    parser.add_argument("--exact-recovery", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--localization", type=Path, required=True)
    parser.add_argument("--raisd", type=Path, required=True)
    parser.add_argument("--figure", type=Path, action="append", default=[])
    parser.add_argument("--job-id", action="append", default=[])
    parser.add_argument("--final-sealed-sha256", required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--bed-out", type=Path, required=True)
    args = parser.parse_args()

    paths = {
        "generic_recovery": args.generic_recovery,
        "chromosome_recovery": args.chromosome_recovery,
        "exact_recovery": args.exact_recovery,
        "candidate": args.candidate,
        "audit": args.audit,
        "localization": args.localization,
        "raisd": args.raisd,
    }
    loaded = {name: load(path) for name, path in paths.items()}
    gates = validate(
        loaded["generic_recovery"],
        loaded["chromosome_recovery"],
        loaded["exact_recovery"],
        loaded["candidate"],
        loaded["raisd"],
    )
    audit = loaded["audit"]
    candidate = loaded["candidate"]
    start, end = candidate["interval_selection"]["interval"]
    result = {
        "challenge": "Single hard sweep",
        "status": "method frozen; no submission",
        "submission_count": 0,
        "testing_prediction": {
            "chromosome": str(candidate["chromosome"]),
            "start": start,
            "end": end,
            "center": candidate["interval_selection"]["center"],
        },
        "selected_method": {
            "localizer": candidate["selected_method"],
            "interval_policy": "challenge-score calibration plus independent-consensus guard",
            "calibration_half_width": candidate["interval_selection"]["selection_optimal_half_width"],
            "frozen_half_width": candidate["interval_selection"]["final_half_width"],
        },
        "validation": gates,
        "observed_estimates": candidate["observed_estimates"],
        "data_provenance": {
            "testing_vcf_sha256": audit["files"]["vcf"]["sha256"],
            "accessibility_bed_sha256": audit["files"]["accessible_bed"]["sha256"],
            "final_vcf_sealed_sha256": args.final_sealed_sha256,
            "accepted_biallelic_sites": audit["accepted_biallelic_sites"],
            "haplotypes": audit["haplotypes"],
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
    args.bed_out.parent.mkdir(parents=True, exist_ok=True)
    args.bed_out.write_text(f"{candidate['chromosome']}\t{start}\t{end}\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
