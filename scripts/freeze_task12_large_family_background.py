from __future__ import annotations

import csv
import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

from ghist_relatedness.labels import HUMAN_LABELS
from ghist_submission.validate import all_pairs, vcf_samples, validate_submission


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT_ROOT = ROOT / "artifacts/task12_large_family_background"
DEVELOPMENT = ARTIFACT_ROOT / "empirical_background_v1/summary.json"
HOLDOUT = ARTIFACT_ROOT / "empirical_background_v1/holdout_summary.json"
CONSENSUS = ARTIFACT_ROOT / "empirical_background_v1_guarded/summary.json"
CANDIDATE = ARTIFACT_ROOT / "empirical_background_v1_guarded/candidate.csv"
FREEZE = ROOT / "experiments/large_family_background/task12_freeze.json"
FIGURE_ROOT = ROOT / "figures/task12_large_family_background"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def relative(path: Path) -> str:
    return str(path.relative_to(ROOT))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-vcf", type=Path, required=True)
    args = parser.parse_args()
    cohort = vcf_samples(args.sample_vcf)
    coverage = validate_submission(CANDIDATE, 12, "testing", all_pairs(cohort))
    previous = json.loads(FREEZE.read_text())
    development = json.loads(DEVELOPMENT.read_text())
    holdout = json.loads(HOLDOUT.read_text())
    consensus = json.loads(CONSENSUS.read_text())
    with CANDIDATE.open() as handle:
        rows = list(csv.DictReader(handle))
    accepted = set(HUMAN_LABELS.values())
    pairs = [tuple(sorted((row["indID1"], row["indID2"]))) for row in rows]
    samples = {sample for pair in pairs for sample in pair}
    labels = Counter(row["relation"] for row in rows)
    checks = {
        "empirical_background_holdout_gate": bool(holdout["gate"]["passed"]),
        "all_input_pairs_present": coverage["input_pair_coverage_verified"],
        "pairs_unique": len(set(pairs)) == len(pairs),
        "all_input_individuals_present": samples == set(cohort),
        "all_labels_accepted": set(labels).issubset(accepted),
        "official_relation_only_schema": set(rows[0])
        == {"indID1", "indID2", "relation"},
        "topology_cannot_override_meiosis_degree": bool(
            consensus["require_flat_degree_agreement"]
        ),
        "submission_count_zero": consensus["submission_count"] == 0,
    }
    checks["passed"] = all(checks.values())
    if not checks["passed"]:
        raise RuntimeError(f"Task 12 freeze gate failed: {checks}")

    source_paths = (DEVELOPMENT, HOLDOUT, CONSENSUS, CANDIDATE)
    payload = {
        "challenge": "Relatedness in a large human family with background relatedness",
        "status": "FROZEN",
        "betty_job_ids": ["7828946", "7828972", "7829025"],
        "data_provenance": {
            **previous["data_provenance"],
            "sample_vcf_sha256": sha256(args.sample_vcf),
        },
        "figures": [
            relative(FIGURE_ROOT / "task12_validation.pdf"),
            relative(FIGURE_ROOT / "task12_validation.png"),
        ],
        "selected_method": {
            "pair_model": "5 cM histogram-gradient classifier using corrected segment counts, total sharing, and length-bin features",
            "background_control": "background length-bin rates estimated from unlabeled Task 12 versus Task 10 cohort distributions",
            "selection_rule": "selected on empirical short-segment contamination regimes and promoted through 20 untouched paired comparisons",
            "topology_consensus": "branch reconstruction may refine relation type only when its meiosis degree agrees with the calibrated pair model",
        },
        "source_artifacts": {relative(path): sha256(path) for path in source_paths},
        "validation": {
            "input_pair_coverage": coverage,
            "development_empirical_background_exact_accuracy": development[
                "development_evaluation"
            ]["mean_exact_accuracy"],
            "development_empirical_background_meioses_accuracy": development[
                "development_evaluation"
            ]["mean_meioses_accuracy"],
            "untouched_empirical_background_replicates": holdout["holdout_replicates"],
            "untouched_empirical_background_comparisons": holdout["comparisons"],
            "previous_exact_accuracy_on_untouched_empirical_background": holdout[
                "baseline"
            ]["mean_exact_accuracy"],
            "frozen_exact_accuracy_on_untouched_empirical_background": holdout[
                "candidate"
            ]["mean_exact_accuracy"],
            "mean_paired_exact_gain": holdout["mean_paired_exact_gain"],
            "minimum_paired_exact_gain": holdout["minimum_paired_exact_gain"],
            "positive_exact_comparisons": holdout["positive_exact_comparisons"],
            "previous_meioses_accuracy_on_untouched_empirical_background": holdout[
                "baseline"
            ]["mean_meioses_accuracy"],
            "frozen_meioses_accuracy_on_untouched_empirical_background": holdout[
                "candidate"
            ]["mean_meioses_accuracy"],
            "mean_paired_meioses_gain": holdout["mean_paired_meioses_gain"],
            "degree_guarded_topology_refinements": consensus["source_counts"].get(
                "background-guarded accepted pedigree tuple", 0
            ),
            "robustness_gate": checks,
        },
        "numerical_diagnostics": {
            "rejected_background_augmented_raw_model": "failed development comparison",
            "rejected_pairwise_segment_marginalization": "improved light contamination but regressed clean and heavy regimes",
            "rejected_long_only_validation": "discarded after unlabeled cohort audit showed dominant 5-20 cM background",
            "empirical_background_profile_by_length_bin": development[
                "background_profile"
            ],
        },
        "testing_prediction": {
            "pairs": len(rows),
            "individuals": len(samples),
            "path": relative(CANDIDATE),
            "schema": ["indID1", "indID2", "relation"],
            "label_counts": dict(labels),
        },
        "submission_count": 0,
    }
    FREEZE.write_text(json.dumps(payload, indent=2) + "\n")
    print(
        json.dumps(
            {"status": payload["status"], "validation": payload["validation"]}, indent=2
        )
    )


if __name__ == "__main__":
    main()
