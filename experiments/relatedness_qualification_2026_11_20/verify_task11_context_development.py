"""Independently score the retained fixed family-context development experiment."""
import hashlib
import json
import math
import tarfile
from datetime import UTC, datetime
from pathlib import Path

from scripts.score_relatedness_cohorts import evaluate_cohorts
from scripts.score_relatedness_published_metric import evaluate_published_metric


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "experiments/task11_development_2026_11_20"
OUT = ROOT / "experiments/relatedness_qualification_2026_11_20/task11_context_development"
ARCHIVE_SHA = "764e343977da10328a01f5036d2454cf7abfd940bec6c2012ad326d4f30088e7"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    archive_path = BASE / "context_diagnostic_v1_compact.tar.gz"
    assert digest(archive_path) == ARCHIVE_SHA
    folder = BASE / "context_diagnostic_v1"
    before = {p.name: digest(p) for p in folder.iterdir() if p.is_file()}
    with tarfile.open(archive_path) as archive:
        members = archive.getmembers()
        assert len(members) == len({m.name for m in members})
        for member in members:
            path = Path(member.name)
            assert not path.is_absolute() and ".." not in path.parts
            assert not member.issym() and not member.islnk()
            if member.isfile():
                assert digest(BASE / path) == hashlib.sha256(archive.extractfile(member).read()).hexdigest()
    report = json.loads((folder / "report.json").read_text())
    for name, actual in before.items():
        if name != "report.json":
            assert report["outputs_sha256"][name] == actual, name
    assert digest(folder / "prospective_plan.json") == "40ba7775182af3e3f1c05a14be2db84793b30b507d90592f4c848be61fd905a3"
    results = {}
    for split in ("calibration", "exposed"):
        raw = json.loads((folder / f"scoring_{split}.json").read_text())
        assert len(raw["cohorts"]) == (4 if split == "calibration" else 8)
        if split == "exposed":
            for index, cohort in enumerate(raw["cohorts"]):
                old = BASE / "evidence_v3/results_v3" / f"confirmation_{index:02d}/truth.json"
                assert cohort["truth"] == json.loads(old.read_text())["pairs"]
        scores = evaluate_cohorts(raw["cohorts"], raw["labels"], bootstrap_seed=20301120)
        metric = evaluate_published_metric(raw["cohorts"], raw["labels"])
        claimed = report["reports"][split]["context"]
        assert scores["pooled"]["pairs"] == claimed["truth_pairs"]
        assert scores["pooled"]["correct"] == claimed["correct"]
        assert scores["missing_pairs"] == claimed["missing_true_pairs"]
        assert scores["extra_pairs"] == claimed["extra_cross_family_pairs"]
        for label, row in scores["pooled"]["per_label"].items():
            assert row["pairs"] == claimed["counts"][label]
            assert row["correct"] == claimed["correct_counts"][label]
        for actual, expected in ((scores["pooled"]["exact_accuracy"], claimed["accuracy"]),
                                 (scores["pooled"]["balanced_accuracy"], claimed["balanced_accuracy"]),
                                 (scores["calibration"]["ece"], claimed["ece_10_equal_width"])):
            assert math.isclose(actual, expected, abs_tol=1e-12, rel_tol=0)
        results[split] = {"scores": scores, "published_metric": metric}
    assert before == {p.name: digest(p) for p in folder.iterdir() if p.is_file()}
    OUT.mkdir(exist_ok=False)
    summary = {}
    for split, result in results.items():
        for name, data in result.items():
            (OUT / f"{split}_{name}.json").write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
        s = result["scores"]
        summary[split] = {"correct": s["pooled"]["correct"], "pairs": s["pooled"]["pairs"],
                          "accuracy": s["pooled"]["exact_accuracy"], "balanced_accuracy": s["pooled"]["balanced_accuracy"],
                          "minimum_label_recall": min(row["recall"] for row in s["pooled"]["per_label"].values()),
                          "weighted_jaccard": result["published_metric"]["pooled"]["weighted_jaccard"],
                          "ece": s["calibration"]["ece"], "bootstrap_lower": s["bootstrap"]["balanced_accuracy_interval"][0]}
    audit = {"verified_at_utc": datetime.now(UTC).isoformat(), "archive_sha256": ARCHIVE_SHA,
             "input_sha256_before_and_after": before, "audit_source_sha256": digest(Path(__file__)),
             "summary": summary, "scope": "Independent raw development rescoring; calibration/exposed cohorts are not fresh confirmation.",
             "release_authorized": False}
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2, allow_nan=False) + "\n")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
