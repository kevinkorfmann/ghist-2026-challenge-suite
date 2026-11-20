"""Independently check bound Task 10 development predictions and rescore them."""
import argparse
import hashlib
import itertools
import json
from datetime import UTC, datetime
from pathlib import Path

from scripts.score_relatedness_cohorts import evaluate_cohorts
from scripts.score_relatedness_published_metric import evaluate_published_metric


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def pair(row):
    return tuple(sorted((row["indID1"], row["indID2"])))


def main(args):
    results, bindings, remote_bindings, manifest_hashes = {}, {}, {}, {}
    for relative in args.manifest:
        manifest_path = args.base / relative
        manifest_hashes[relative] = sha(manifest_path)
        manifest = json.loads(manifest_path.read_text())
        assert "DEVELOPMENT" in manifest["stage"]
        for name, expected in manifest["files_sha256"].items():
            path = args.base / name
            assert not Path(name).is_absolute() and ".." not in Path(name).parts
            assert sha(path) == expected, name
            assert bindings.get(name, expected) == expected
            bindings[name] = expected
        for name, row in manifest.get("remote_bound_files", {}).items():
            assert sha(args.base / name) == row["sha256"], name
            assert remote_bindings.get(row["remote_path"], row["sha256"]) == row["sha256"]
            remote_bindings[row["remote_path"]] = row["sha256"]
        raw = json.loads((args.base / manifest["scorer_input"]).read_text())
        assert len(raw["labels"]) == len(set(raw["labels"])) == 14
        assert len(raw["cohorts"]) == 8
        prediction_paths = [name for name in manifest["files_sha256"]
                            if name.endswith("predictions.json") and "/cohort_" not in name]
        if prediction_paths:
            assert len(prediction_paths) == 1
            predictions = json.loads((args.base / prediction_paths[0]).read_text())
        else:
            assert manifest["variant"] == "baseline"
            predictions = None
        for cohort in raw["cohorts"]:
            index = int(cohort["cohort_id"].rsplit("_", 1)[1])
            assert 40 <= index <= 47
            prefix = args.base / f"detected_panel_evidence/cohort_{index:03d}"
            truth = json.loads((prefix / "truth.json").read_text())
            identifiers = json.loads((prefix / "identifiers.json").read_text())
            assert len(set(identifiers.values())) == len(identifiers)
            expected_pairs = set(itertools.combinations(sorted(identifiers.values()), 2))
            assert truth == cohort["truth"]
            assert {pair(row) for row in truth} == expected_pairs
            assert {pair(row) for row in cohort["predictions"]} == expected_pairs
            assert len(cohort["predictions"]) == len(truth) == len(expected_pairs)
            if predictions is None:
                bound_predictions = json.loads((prefix / "baseline_predictions.json").read_text())
            else:
                bound_predictions = [{k: v for k, v in row.items() if k != "cohort_id"}
                                     for row in predictions if row["cohort_id"] == cohort["cohort_id"]]
            assert cohort["predictions"] == bound_predictions
        scores = evaluate_cohorts(raw["cohorts"], raw["labels"], bootstrap_seed=20301120)
        published = evaluate_published_metric(raw["cohorts"], raw["labels"])
        assert scores["pooled"]["pairs"] == 585
        assert scores["missing_pairs"] == scores["extra_pairs"] == 0
        key = Path(relative).parent.name + "_" + manifest["variant"]
        assert key not in results
        pooled = scores["pooled"]
        results[key] = {"scores": scores, "published_metric": published,
                        "summary": {"correct": pooled["correct"], "pairs": pooled["pairs"],
                                    "accuracy": pooled["exact_accuracy"],
                                    "macro_recall": pooled["balanced_accuracy"],
                                    "minimum_label_recall": min(row["recall"] for row in pooled["per_label"].values()),
                                    "weighted_jaccard": published["pooled"]["weighted_jaccard"],
                                    "ece": scores["calibration"]["ece"] if scores["calibration"] else None,
                                    "bootstrap": scores["bootstrap"]}}
    for name, expected in bindings.items():
        assert sha(args.base / name) == expected, name
    for name, expected in manifest_hashes.items():
        assert sha(args.base / name) == expected, name
    args.output.mkdir(exist_ok=False)
    for key, result in results.items():
        for metric in ("scores", "published_metric"):
            (args.output / f"{key}_{metric}.json").write_text(json.dumps(result[metric], indent=2, allow_nan=False) + "\n")
    report = {"verified_at_utc": datetime.now(UTC).isoformat(), "manifest_sha256": manifest_hashes,
              "bound_file_count": len(bindings), "files_sha256": bindings,
              "source_sha256": sha(Path(__file__)), "remote_bindings_to_verify": remote_bindings,
              "scope": "Independent raw-row, complete-roster, pair-coverage and before/after local hash audit. Remote bindings require the separate remote report. Development only.",
              "summary": {key: result["summary"] for key, result in results.items()},
              "release_authorized": False}
    (args.output / "audit.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps(report["summary"]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--manifest", action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    main(parser.parse_args())
