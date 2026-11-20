"""Freeze only the new full-input Task 7 method with complete matched recovery.

Archived caller-region/unmasked artifacts are preserved but are not accepted.
No matched biological recovery producer exists yet; current candidates remain
provisional. A future recovery must satisfy the artifact contract below.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

from ghist_single_sweep.background import METHOD, SETTINGS, predict_masked_ten
from ghist_single_sweep.peaks import read_raisd
from ghist_single_sweep.scan import read_bed_intervals, read_phased_biallelic_vcf
from ghist_submission.validate import validate_submission
from scripts.predict_task7_sweep_bgs import (
    sha256,
    source_manifest,
    validate_raisd_site_count,
)

SCHEMA = "task7_matched_complete_recovery_v1"
ARTIFACTS = {
    "masked.vcf",
    "accessible.bed",
    "RAiSD_Report.task7",
    "RAiSD_Info.task7",
    "candidate.bed",
}


def verify_candidate(path):
    """Recompute the exact new predictor from retained observation artifacts."""
    payload = json.loads(path.read_text())
    if payload.get("method") != METHOD or payload.get("settings") != SETTINGS:
        raise ValueError("Legacy or mismatched Task 7 method/settings")
    if payload.get("source_sha256") != source_manifest():
        raise ValueError("Stale or mismatched Task 7 source provenance")
    hashes = payload.get("artifact_sha256", {})
    if set(hashes) != ARTIFACTS or any(
        sha256(path.parent / name) != digest for name, digest in hashes.items()
    ):
        raise ValueError("Missing or mismatched Task 7 observation artifacts")
    provenance = payload.get("input_provenance", {})
    if set(provenance) != {"vcf_sha256", "accessibility_sha256"} or any(
        not isinstance(v, str) or re.fullmatch(r"[0-9a-f]{64}", v) is None
        for v in provenance.values()
    ):
        raise ValueError("Missing Task 7 original input provenance")
    dataset = payload.get("dataset")
    if dataset not in ("testing", "final"):
        raise ValueError("Task 7 dataset identity missing")
    expected = "7" if dataset == "testing" else "9"
    variants = read_phased_biallelic_vcf(path.parent / "masked.vcf")
    if variants.contig.removeprefix("chr") != expected:
        raise ValueError("Task 7 VCF chromosome differs from dataset")
    validate_raisd_site_count(path.parent / "RAiSD_Info.task7", len(variants.positions))
    mask = read_bed_intervals(path.parent / "accessible.bed", variants.contig)
    rebuilt = predict_masked_ten(
        variants, mask, read_raisd(path.parent / "RAiSD_Report.task7")
    )
    if any(payload.get(key) != value for key, value in rebuilt.items()):
        raise ValueError("Task 7 candidate differs from recomputed observations")
    lo, hi = rebuilt["interval"]
    if (path.parent / "candidate.bed").read_text() != f"{expected}\t{lo}\t{hi}\n":
        raise ValueError("Task 7 BED differs from recomputed interval")
    validate_submission(path.parent / "candidate.bed", 7, dataset)
    return payload


def validate_matched_recovery(recovery):
    if recovery.get("schema") != SCHEMA:
        raise ValueError("Legacy Task 7 recovery is not new matched complete recovery")
    if (
        recovery.get("status") != "complete_independent_confirmation"
        or recovery.get("public_testing_answers_used") is not False
    ):
        raise ValueError(
            "Task 7 requires new independent confirmation, not development/reanalysis"
        )


def score_recovery(recovery_path, candidate):
    recovery = json.loads(recovery_path.read_text())
    validate_matched_recovery(recovery)
    design_path = recovery_path.parent / recovery["design_path"]
    design_hash = sha256(design_path)
    if design_hash != recovery.get("design_sha256"):
        raise ValueError("Mismatched prespecified recovery design")
    design = json.loads(design_path.read_text())
    if (
        design.get("status") != "prespecified"
        or design.get("method") != METHOD
        or design.get("settings") != SETTINGS
        or design.get("source_sha256") != source_manifest()
    ):
        raise ValueError(
            "Recovery design has stale method/source or is not prespecified"
        )
    context = design.get("simulation_context", {})
    if (
        context.get("background_selection") is not True
        or context.get("slim_scaling_factor") != 1
        or not context.get("demography")
        or not context.get("dfe")
        or not context.get("simulation_source_sha256")
    ):
        raise ValueError(
            "Recovery requires explicit exact-scale BGS simulation provenance"
        )
    exposure = {
        k: candidate[k]
        for k in ("chromosome", "contig_length", "haplotypes", "callable_bases")
    }
    exposure["canonical_mask_sha256"] = candidate["artifact_sha256"]["accessible.bed"]
    if design.get("exposure") != exposure or design.get("software") != candidate.get(
        "software"
    ):
        raise ValueError("Recovery exposure/software does not match candidate")
    groups = {
        part: design.get(f"{part}_seeds", []) for part in ("calibration", "evaluation")
    }
    seeds = groups["calibration"] + groups["evaluation"]
    if (
        any(not group for group in groups.values())
        or any(type(seed) is not int or seed < 0 for seed in seeds)
        or len(set(seeds)) != len(seeds)
    ):
        raise ValueError(
            "Calibration/evaluation need nonempty disjoint unique seed groups"
        )
    records = recovery.get("records", [])
    if len(records) != len(seeds) or {r.get("seed") for r in records} != set(seeds):
        raise ValueError("Recovery must contain every planned seed exactly once")
    evaluations = []
    tree_hashes = set()
    input_hashes = {candidate["artifact_sha256"]["masked.vcf"]}
    for record in records:
        seed = record["seed"]
        part = "calibration" if seed in groups["calibration"] else "evaluation"
        if record.get("partition") != part or record.get("state") != "COMPLETED":
            raise ValueError("Recovery partitions/state differ from complete design")
        prediction_path = recovery_path.parent / record["prediction_path"]
        truth_path = recovery_path.parent / record["truth_path"]
        if sha256(prediction_path) != record.get("prediction_sha256") or sha256(
            truth_path
        ) != record.get("truth_sha256"):
            raise ValueError("Recovery record source hash mismatch")
        prediction = verify_candidate(prediction_path)
        input_hash = prediction["artifact_sha256"]["masked.vcf"]
        if input_hash in input_hashes:
            raise ValueError("Recovery reuses candidate or another replicate VCF")
        input_hashes.add(input_hash)
        if (
            any(
                prediction[k] != exposure[k]
                for k in ("chromosome", "contig_length", "haplotypes", "callable_bases")
            )
            or prediction["artifact_sha256"]["accessible.bed"]
            != exposure["canonical_mask_sha256"]
            or prediction.get("software") != candidate.get("software")
        ):
            raise ValueError("Recovery prediction exposure/software mismatch")
        truth = json.loads(truth_path.read_text())
        if (
            truth.get("seed") != seed
            or truth.get("design_sha256") != design_hash
            or truth.get("simulation_context") != context
        ):
            raise ValueError("Recovery truth identity/context mismatch")
        position = truth.get("position")
        if type(position) is not int or not 0 <= position < prediction["contig_length"]:
            raise ValueError("Recovery truth coordinate invalid")
        tree_path = truth_path.parent / truth["tree_path"]
        tree_hash = sha256(tree_path)
        if tree_hash != truth.get("tree_sha256"):
            raise ValueError("Recovery simulation tree hash mismatch")
        if tree_hash in tree_hashes:
            raise ValueError("Recovery reuses a simulation tree across seeds")
        tree_hashes.add(tree_hash)
        if part == "evaluation":
            lo, hi = prediction["interval"]
            evaluations.append(
                {
                    "seed": seed,
                    "eligible": prediction["consensus_eligible"],
                    "covered": lo <= position < hi,
                    "absolute_error": abs(prediction["center"] - position),
                }
            )
    requirements = design.get("requirements", {})
    for name in (
        "minimum_evaluation_replicates",
        "minimum_eligible_evaluation_replicates",
    ):
        if type(requirements.get(name)) is not int or requirements[name] < 2:
            raise ValueError("Prespecified replicate requirements must be integers >=2")
    minimum = requirements.get("minimum_conditional_coverage_lower95")
    if (
        not isinstance(minimum, (int, float))
        or not math.isfinite(minimum)
        or not 0 < minimum <= 1
    ):
        raise ValueError("Prespecified coverage requirement must be finite in (0,1]")
    eligible = [r for r in evaluations if r["eligible"]]
    covered = sum(r["covered"] for r in eligible)
    lower = (
        max(
            0.0, covered / len(eligible) - math.sqrt(math.log(40) / (2 * len(eligible)))
        )
        if eligible
        else 0.0
    )
    gates = {
        "evaluation_count": len(evaluations)
        >= requirements["minimum_evaluation_replicates"],
        "eligible_count": len(eligible)
        >= requirements["minimum_eligible_evaluation_replicates"],
        "conditional_coverage_lower95": lower >= minimum,
        "candidate_eligible": candidate["consensus_eligible"],
    }
    if not all(gates.values()):
        raise ValueError(f"Task 7 prespecified recovery requirements fail: {gates}")
    return {
        "gates": gates,
        "evaluation_records": evaluations,
        "conditional_coverage_lower95": lower,
        "coverage_assumption": "Independent evaluation replicates conditional on the fixed eligibility rule; conservative Hoeffding lower bound for average coverage, allowing heterogeneous cases",
        "design_sha256": design_hash,
        "recovery_sha256": sha256(recovery_path),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--matched-recovery", type=Path, required=True)
    parser.add_argument("--bed", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.bed.exists() or args.out.exists():
        raise ValueError("Refusing to overwrite Task 7 freeze artifacts")
    candidate = verify_candidate(args.candidate)
    validation = score_recovery(args.matched_recovery, candidate)
    args.bed.parent.mkdir(parents=True, exist_ok=True)
    args.bed.write_bytes((args.candidate.parent / "candidate.bed").read_bytes())
    validate_submission(args.bed, 7, candidate["dataset"])
    result = {
        "challenge": 7,
        "method": METHOD,
        "status": "method package frozen; scientific readiness not established",
        "scientifically_validated": False,
        "release_ready": False,
        "required_review": "A separately reviewed scientific protocol is required; declared performance requirements are not scientific approval",
        "submission_count": 0,
        "candidate_sha256": sha256(args.candidate),
        "bed_sha256": sha256(args.bed),
        "validation": validation,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
