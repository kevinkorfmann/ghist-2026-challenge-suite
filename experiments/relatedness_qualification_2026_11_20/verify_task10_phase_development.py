"""Independently rescore Task10 phased development and its matched controls."""
import hashlib
from itertools import combinations
import json
from pathlib import Path

from scripts.score_relatedness_cohorts import evaluate_cohorts
from scripts.score_relatedness_published_metric import evaluate_published_metric


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "experiments/task10_development_2026_11_20"
OUT = ROOT / "experiments/relatedness_qualification_2026_11_20/task10_phase_development"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path):
    return json.loads(path.read_text())


def main():
    reports = {}
    for variant in ("baseline", "phase", "phase_free"):
        manifest_path = BASE / "phase_scorer_inputs" / f"{variant}.manifest.json"
        manifest = load(manifest_path)
        for name, expected in manifest["files_sha256"].items():
            assert digest(BASE / name) == expected, name
        data = load(BASE / manifest["scorer_input"])
        for index, cohort in enumerate(data["cohorts"]):
            folder = BASE / "phase_evidence" / f"phase_development_cohort_{index}"
            truth = load(folder / "truth.json")
            if variant == "baseline":
                prediction_path = folder / "representation_comparison/raw_haplotype.json"
            else:
                prediction_path = BASE / "phase_evidence/phase_predictions" / f"cohort_{index}" / variant / "predictions.json"
            assert cohort["truth"] == truth
            assert cohort["predictions"] == load(prediction_path)
            roster = {sample for row in truth for sample in (row["indID1"], row["indID2"])}
            pairs = {tuple(sorted((row["indID1"], row["indID2"]))) for row in truth}
            assert pairs == set(combinations(sorted(roster), 2))
            assert len(pairs) == len(truth)
        score = evaluate_cohorts(data["cohorts"], data["labels"])
        published = evaluate_published_metric(data["cohorts"], data["labels"])
        for name, expected in manifest["files_sha256"].items():
            assert digest(BASE / name) == expected, name
        reports[variant] = {
            "raw_cohort_scores": score, "published_metric": published,
            "verified_manifest_sha256": digest(manifest_path),
            "input_sha256_before_and_after": manifest["files_sha256"],
            "calibration_scope": "Raw model probabilities, before the separately selected development temperature; not confirmation",
            "audit_source_sha256": digest(Path(__file__)),
            "scorer_sha256": digest(ROOT / "scripts/score_relatedness_cohorts.py"),
            "published_metric_scorer_sha256": digest(ROOT / "scripts/score_relatedness_published_metric.py"),
        }
    OUT.mkdir(exist_ok=False)
    summary = {}
    for variant, report in reports.items():
        (OUT / f"{variant}.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
        score = report["raw_cohort_scores"]
        summary[variant] = {"correct": score["pooled"]["correct"], "pairs": score["pooled"]["pairs"],
                            "accuracy": score["pooled"]["exact_accuracy"], "macro": score["pooled"]["balanced_accuracy"],
                            "weighted_jaccard": report["published_metric"]["pooled"]["weighted_jaccard"],
                            "raw_ece": score["calibration"]["ece"], "coverage_passed": score["pair_coverage_passed"]}
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
