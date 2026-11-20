"""Rebuild Task11 v3 scores from retained raw rows, checking their bindings."""
import hashlib
from itertools import combinations
import json
from pathlib import Path

from scripts.score_relatedness_cohorts import evaluate_cohorts
from scripts.score_relatedness_published_metric import evaluate_published_metric


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "experiments/task11_development_2026_11_20"
EVIDENCE = BASE / "evidence_v3"
OUT = ROOT / "experiments/relatedness_qualification_2026_11_20/task11_confirmation_v3"
PROTOCOL_HASH = "f5ec05f6a4e3d27d82021bf6d1623e1392735f1896d9a2bda2249020b7120996"
CANDIDATE_HASH = "c966d2f676057ebaf60b862b73f2b69d076e92a7574457497e339f855c262f50"
ARCHIVE_HASH = "c99b2dfa7fee4447111cb3e43b45bad109f595f05ba204754ccefdc137cfd573"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path):
    return json.loads(path.read_text())


def pairs(rows):
    return {tuple(sorted((row["indID1"], row["indID2"]))) for row in rows}


def main():
    before = {str(path.relative_to(EVIDENCE)): digest(path) for path in EVIDENCE.rglob("*.json")}
    assert digest(BASE / "evidence_v3.tar.gz") == ARCHIVE_HASH
    assert before["protocol_v3.json"] == PROTOCOL_HASH
    assert before["candidate_v3/candidate.json"] == CANDIDATE_HASH
    protocol = load(EVIDENCE / "protocol_v3.json")
    summary = load(EVIDENCE / "summary_v3.json")
    candidate = load(EVIDENCE / "candidate_v3/candidate.json")
    crosscheck = load(EVIDENCE / "independent_scoring_input.json")
    assert summary["protocol_sha256"] == candidate["protocol_sha256"] == PROTOCOL_HASH
    assert summary["candidate_sha256"] == CANDIDATE_HASH
    assert candidate["source_sha256"] == protocol["source_sha256"]
    assert before["independent_scoring_input.json"] == summary["crosscheck_input_sha256"]
    assert crosscheck["labels"] == protocol["labels"]
    expected_configs = [row for row in protocol["cohorts"] if row["split"] == "confirmation"]
    cohorts, partitions = [], []
    for index, config in enumerate(expected_configs):
        folder = EVIDENCE / "results_v3" / f"confirmation_{index:02d}"
        result, cohort_protocol = load(folder / "result.json"), load(folder / "protocol.json")
        remote = next(key for key in summary["confirmation_results_sha256"] if key.endswith(f"confirmation_{index:02d}/result.json"))
        assert digest(folder / "result.json") == summary["confirmation_results_sha256"][remote]
        assert cohort_protocol["master_sha256"] == PROTOCOL_HASH
        assert cohort_protocol["candidate_sha256"] == CANDIDATE_HASH
        assert cohort_protocol["source_sha256"] == protocol["source_sha256"]
        assert cohort_protocol["cohort"] == result["cohort"] == config
        assert digest(folder / "protocol.json") == result["protocol_sha256"]
        for path in folder.glob("*.json"):
            if path.name != "result.json":
                assert digest(path) == result["output_sha256"][path.name], path
        truth, predictions, clusters = load(folder / "truth.json"), load(folder / "predictions.json"), load(folder / "clusters.json")
        actual_families = {frozenset(family) for family in truth["families"]}
        inferred_families = {frozenset(family) for family in clusters["families"]}
        expected_pairs = {tuple(sorted(pair)) for family in actual_families for pair in combinations(family, 2)}
        assert pairs(truth["pairs"]) == expected_pairs
        assert len(truth["pairs"]) == len(expected_pairs)
        roster = truth["samples"]
        for partition in (truth["families"], clusters["families"]):
            flattened = [sample for family in partition for sample in family]
            assert len(flattened) == len(set(flattened)) == len(roster)
            assert set(flattened) == set(roster)
        row = {"cohort_id": str(config["seed"]), "truth": truth["pairs"], "predictions": predictions}
        assert row == crosscheck["cohorts"][index]
        cohorts.append(row)
        partitions.append({"cohort_id": row["cohort_id"], "truth_families": len(actual_families),
                           "exact_families": len(actual_families & inferred_families),
                           "truth_pairs": len(expected_pairs), "predicted_pairs": len(predictions),
                           "exact_partition": actual_families == inferred_families})
    report = evaluate_cohorts(cohorts, protocol["labels"], bootstrap_seed=20301120)
    published = evaluate_published_metric(cohorts, protocol["labels"])
    after = {str(path.relative_to(EVIDENCE)): digest(path) for path in EVIDENCE.rglob("*.json")}
    assert before == after
    gates = protocol["gates"]
    checks = {
        "accuracy": report["pooled"]["exact_accuracy"] >= gates["accuracy_min"],
        "balanced_accuracy": report["pooled"]["balanced_accuracy"] >= gates["balanced_accuracy_min"],
        "each_label_recall": min(row["recall"] for row in report["pooled"]["per_label"].values()) >= gates["each_label_recall_min"],
        "bootstrap_lower": report["bootstrap"]["balanced_accuracy_interval"][0] >= gates["bootstrap_balanced_lower95_min"],
        "family_recovery": sum(row["exact_families"] for row in partitions) / sum(row["truth_families"] for row in partitions) >= gates["exact_family_recovery_min"],
        "missing_samples": True,
        "missing_true_pairs": report["missing_pairs"] <= gates["missing_true_pairs_max"],
        "extra_pairs": report["extra_pairs"] <= gates["extra_cross_family_pairs_max"],
        "calibration": report["calibration"]["ece"] <= gates["expected_calibration_error_max"],
    }
    assert checks == summary["checks"]
    OUT.mkdir(exist_ok=False)
    evidence = {"protocol_sha256": PROTOCOL_HASH, "candidate_sha256": CANDIDATE_HASH,
                "archive_sha256": ARCHIVE_HASH, "input_sha256_before_and_after": before,
                "partitions": partitions, "independent_checks": checks,
                "predeclared_criteria_passed": all(checks.values()), "release_authorized": False,
                "audit_source_sha256": digest(Path(__file__)),
                "scorer_sha256": digest(ROOT / "scripts/score_relatedness_cohorts.py"),
                "published_metric_scorer_sha256": digest(ROOT / "scripts/score_relatedness_published_metric.py")}
    for name, value in [("cohort_scores.json", report), ("published_metric.json", published), ("audit.json", evidence)]:
        (OUT / name).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"correct": report["pooled"]["correct"], "pairs": report["pooled"]["pairs"],
                      "accuracy": report["pooled"]["exact_accuracy"], "balanced_accuracy": report["pooled"]["balanced_accuracy"],
                      "bootstrap": report["bootstrap"]["balanced_accuracy_interval"], "ece": report["calibration"]["ece"],
                      "weighted_jaccard": published["pooled"]["weighted_jaccard"], "checks": checks}))


if __name__ == "__main__":
    main()
