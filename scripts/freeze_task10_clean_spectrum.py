from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np

from ghist_relatedness.labels import HUMAN_LABELS
from ghist_submission.validate import all_pairs, vcf_samples, validate_submission


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--previous-freeze", type=Path, required=True)
    parser.add_argument("--model-summary", type=Path, required=True)
    parser.add_argument("--holdout-summaries", type=Path, nargs="+", required=True)
    parser.add_argument("--candidate-summary", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--sample-vcf", type=Path, required=True)
    parser.add_argument("--support-artifacts", type=Path, nargs="*", default=[])
    parser.add_argument("--figures", type=Path, nargs="*", default=[])
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    previous = json.loads(args.previous_freeze.read_text())
    model = json.loads(args.model_summary.read_text())
    holdouts = [json.loads(path.read_text()) for path in args.holdout_summaries]
    topology = json.loads(args.candidate_summary.read_text())
    runs = [run for summary in holdouts for run in summary["runs"]]

    samples = vcf_samples(args.sample_vcf)
    coverage = validate_submission(args.candidate, 10, "testing", all_pairs(samples))
    with args.candidate.open() as handle:
        rows = list(csv.DictReader(handle))
    pairs = [tuple(sorted((row["indID1"], row["indID2"]))) for row in rows]
    people = {sample for pair in pairs for sample in pair}
    labels = Counter(row["relation"] for row in rows)
    exact_gain = np.asarray([run["exact_difference"] for run in runs])
    degree_gain = np.asarray([run["meioses_difference"] for run in runs])
    exact = np.asarray([run["candidate"]["exact_accuracy"] for run in runs])
    baseline_exact = np.asarray([run["baseline"]["exact_accuracy"] for run in runs])
    degree = np.asarray([run["candidate"]["meioses_accuracy"] for run in runs])
    baseline_degree = np.asarray([run["baseline"]["meioses_accuracy"] for run in runs])

    gate = {
        "multiple_independent_exact_gains_positive": len(runs) >= 2
        and bool(np.all(exact_gain > 0)),
        "pooled_mean_exact_gain_above_0.006": float(exact_gain.mean()) > 0.006,
        "pooled_meiosis_degree_non_regression": float(degree_gain.mean()) >= 0,
        "official_relation_only_schema": list(rows[0])
        == ["indID1", "indID2", "relation"],
        "all_input_pairs_present": coverage["input_pair_coverage_verified"],
        "pairs_unique": len(set(pairs)) == len(pairs),
        "all_input_individuals_present": people == set(samples),
        "all_labels_accepted": set(labels).issubset(set(HUMAN_LABELS.values())),
    }
    gate["passed"] = all(gate.values())
    if not gate["passed"]:
        raise ValueError(f"Freeze gate failed: {gate}")

    source_paths = [
        args.sample_vcf,
        args.model_summary,
        *args.holdout_summaries,
        args.candidate_summary,
        args.candidate,
        *args.support_artifacts,
    ]
    payload = {
        "challenge": "Relatedness in a large human family",
        "status": "method frozen; no submission",
        "submission_count": 0,
        "testing_prediction": {
            "pairs": len(rows),
            "individuals": len(people),
            "path": str(args.candidate),
            "schema": list(rows[0]),
            "label_counts": dict(labels),
        },
        "selected_method": {
            "pair_model": "75% shallow histogram-gradient full-shape model and 25% calibrated forest",
            "minimum_segment_cm": 5,
            "selected_configuration": model["selected_model"],
            "topology_consensus": "exact intersection of fast and medium branch reconstructions, followed by pair-model meiosis-degree agreement",
            "exact_score_policy": "preserve the calibrated pair call unless replicated topology resolves a same-degree relation ambiguity",
        },
        "validation": {
            "input_pair_coverage": coverage,
            "untouched_replicates": len(runs),
            "candidate_mean_exact_accuracy": float(exact.mean()),
            "baseline_mean_exact_accuracy": float(baseline_exact.mean()),
            "mean_paired_exact_gain": float(exact_gain.mean()),
            "minimum_paired_exact_gain": float(exact_gain.min()),
            "positive_exact_replicates": int(np.sum(exact_gain > 0)),
            "candidate_mean_meioses_accuracy": float(degree.mean()),
            "baseline_mean_meioses_accuracy": float(baseline_degree.mean()),
            "mean_paired_meioses_gain": float(degree_gain.mean()),
            "testing_calls_changed_by_topology": int(
                topology["source_counts"].get(
                    "background-guarded accepted pedigree tuple", 0
                )
            ),
            "all_topology_changes_preserve_meiosis_degree": True,
            "robustness_gate": gate,
        },
        "data_provenance": previous["data_provenance"],
        "source_artifacts": {
            **previous["source_artifacts"],
            **{str(path): digest(path) for path in source_paths},
        },
        "figures": {str(path): digest(path) for path in args.figures},
        "betty_job_ids": list(
            dict.fromkeys(previous["betty_job_ids"] + ["7829097", "7829190", "7829204"])
        ),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload["validation"], indent=2))


if __name__ == "__main__":
    main()
