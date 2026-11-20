"""Independently rescore the failed expanded Task 11 fit and matched ablation."""
import hashlib
import json
import math
import tarfile
from datetime import UTC, datetime
from pathlib import Path

from scripts.score_relatedness_cohorts import evaluate_cohorts
from scripts.score_relatedness_published_metric import evaluate_published_metric


ROOT = Path(__file__).resolve().parents[2]
ARCHIVE = ROOT / "experiments/task11_development_2026_11_20/expanded_context_development_v2_compact.tar.gz"
ARCHIVE_SHA = "00198dbd33b0a8764a2d4656e422295eb28cea992fbd51a9ec9c51c39d014670"
OUT = ROOT / "experiments/relatedness_qualification_2026_11_20/task11_expanded_development"


def main():
    assert hashlib.sha256(ARCHIVE.read_bytes()).hexdigest() == ARCHIVE_SHA
    with tarfile.open(ARCHIVE) as archive:
        members = archive.getmembers()
        assert len(members) == len({m.name for m in members})
        assert all(not Path(m.name).is_absolute() and ".." not in Path(m.name).parts
                   and not m.issym() and not m.islnk() for m in members)
        files = {m.name: archive.extractfile(m).read() for m in members if m.isfile()}

    def load(name):
        return json.loads(files[name])

    def sha(name):
        return hashlib.sha256(files[name]).hexdigest()

    prefix = "expanded_development_audit_v1/"
    fit_prefix = "expanded_context_fit_v2/"
    report = load(prefix + "report.json")
    claimed_audit = load(prefix + "audit.json")
    assert sha(prefix + "report.json") == claimed_audit["fit_report_sha256"]
    for name, expected in claimed_audit["output_sha256"].items():
        assert sha(prefix + name) == expected, name
    exported_fit_bindings = 0
    for name, expected in report["outputs_sha256"].items():
        if fit_prefix + name in files:
            assert sha(fit_prefix + name) == expected, name
            exported_fit_bindings += 1
    assert sha(fit_prefix + "prospective_plan.json") == "0d1c89dc3c6f8925c0c5d88653cf555e28c44343d2f0c98509f9e08fb5f98da8"
    assert sha(prefix + "expanded_training_protocol_v1.json") == "6a775296a270e9808f3d7e906a73d64951d732297c35e56f3f77a496247277e1"
    results, summary = {}, {}
    for split, n in (("calibration", 4), ("exposed", 8)):
        for model in ("base", "context"):
            key = split + "_" + model
            raw = load(prefix + "scoring_" + key + ".json")
            assert len(raw["cohorts"]) == n
            families = 0
            for index, cohort in enumerate(raw["cohorts"]):
                stem = f"{split}_{index:02d}"
                truth = load(prefix + stem + "_truth.json")
                clusters = load(prefix + stem + "_clusters.json")
                assert cohort["truth"] == truth["pairs"]
                truth_partition = {frozenset(group) for group in truth["families"]}
                predicted_partition = {frozenset(group) for group in clusters["families"]}
                assert len(truth_partition) == 8 and predicted_partition == truth_partition
                assert sum(map(len, clusters["families"])) == len(set(truth["samples"]))
                families += 8
            scores = evaluate_cohorts(raw["cohorts"], raw["labels"], bootstrap_seed=20301120)
            metric = evaluate_published_metric(raw["cohorts"], raw["labels"])
            claimed = report["reports"][split][model]
            pooled = scores["pooled"]
            assert pooled["pairs"] == claimed["truth_pairs"]
            assert pooled["correct"] == claimed["correct"]
            assert scores["missing_pairs"] == claimed["missing_true_pairs"] == 0
            assert scores["extra_pairs"] == claimed["extra_cross_family_pairs"] == 0
            for label, row in pooled["per_label"].items():
                assert row["pairs"] == claimed["counts"][label]
                assert row["correct"] == claimed["correct_counts"][label]
            for actual, expected in ((pooled["exact_accuracy"], claimed["accuracy"]),
                                     (pooled["balanced_accuracy"], claimed["balanced_accuracy"]),
                                     (scores["calibration"]["ece"], claimed["ece_10_equal_width"])):
                assert math.isclose(actual, expected, abs_tol=1e-12, rel_tol=0)
            results[key] = {"scores": scores, "published_metric": metric}
            summary[key] = {"correct": pooled["correct"], "pairs": pooled["pairs"],
                            "accuracy": pooled["exact_accuracy"], "balanced_accuracy": pooled["balanced_accuracy"],
                            "minimum_label_recall": min(row["recall"] for row in pooled["per_label"].values()),
                            "weighted_jaccard": metric["pooled"]["weighted_jaccard"],
                            "ece": scores["calibration"]["ece"],
                            "bootstrap_lower": scores["bootstrap"]["balanced_accuracy_interval"][0],
                            "exact_families": families}
    assert claimed_audit["sparse_stress_launch_allowed"] is False
    assert claimed_audit["reserved_confirmation_used"] is False
    OUT.mkdir(exist_ok=False)
    for key, value in results.items():
        for name, data in value.items():
            (OUT / f"{key}_{name}.json").write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    audit = {"verified_at_utc": datetime.now(UTC).isoformat(), "archive_sha256": ARCHIVE_SHA,
             "archive_files": len(files), "exported_fit_bindings_verified": exported_fit_bindings,
             "audit_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             "summary": summary, "scope": "Independent raw-row rescoring and exact family-partition checks on existing development/calibration. Models and omitted raw inputs require separate remote hash audit.",
             "release_authorized": False}
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2, allow_nan=False) + "\n")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
