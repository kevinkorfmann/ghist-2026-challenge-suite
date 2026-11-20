"""Exercise the real Task 9 predictor and submission CLI on fixed simulations."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import subprocess
import sys

from ghist_relatedness.ibd import read_segment_rows, relationship_from_sample
from ghist_relatedness.labels import IBD_LABEL_MEIOSES
from ghist_relatedness.task9_release import sha256


def write_csv(path, columns, rows, encoding="utf-8"):
    with path.open("w", newline="", encoding=encoding) as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--segments", nargs=2, type=Path, required=True)
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--release", type=Path)
    parser.add_argument("--artifact-root", type=Path, required=True)
    parser.add_argument("--sex-map", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    report = {"synthetic_only": True, "public_answers_used": False, "submission_count": 0,
              "model_sha256": sha256(args.model), "datasets": []}
    for index, source in enumerate(args.segments):
        folder = args.out_dir / str(index)
        folder.mkdir(exist_ok=True)
        segments = [row for row in read_segment_rows(source)
                    if relationship_from_sample(row["indID1"]) in IBD_LABEL_MEIOSES]
        # Strip labels from identifiers so the production path cannot recover truth.
        samples = sorted({row[key] for row in segments for key in ("indID1", "indID2")})
        aliases = {sample: f"person_{index}_{number}" for number, sample in enumerate(samples)}
        truth = {}
        for row in segments:
            label = relationship_from_sample(row["indID1"])
            for key in ("indID1", "indID2"):
                row[key] = aliases[row[key]]
            truth[tuple(sorted((row["indID1"], row["indID2"])))] = label
        expected_path = folder / "expected_pairs.csv"
        write_csv(expected_path, ("indID1", "indID2"),
                  [dict(zip(("indID1", "indID2"), pair)) for pair in sorted(truth)])
        baseline_predictions = None
        dataset = {"source_sha256": sha256(source), "expected_pairs": len(truth), "cases": {}}
        for case in ("official_csv", "permuted_reversed_bom_csv"):
            observed = folder / f"{case}.csv"
            columns = list(segments[0])
            rows = segments
            encoding = "utf-8"
            if case != "official_csv":
                rows = [{**row, "indID1": row["indID2"], "indID2": row["indID1"]}
                        for row in reversed(segments)]
                columns = list(reversed(columns))
                encoding = "utf-8-sig"
            write_csv(observed, columns, rows, encoding)
            output = folder / case
            selection = ["--release", str(args.release), "--artifact-root", str(args.artifact_root)] if args.release else ["--model", str(args.model)]
            command = [sys.executable, "scripts/predict_task9_current.py", *selection,
                       "--observed", str(observed), "--sex-map", str(args.sex_map),
                       "--dataset", "testing", "--synthetic-input", "--workers", "2",
                       "--pairs", str(expected_path), "--out-dir", str(output)]
            subprocess.run(command, check=True, stdout=subprocess.DEVNULL)
            candidate = output / "task9_testing_candidate.csv"
            independent_validation = output / "independent_validation.json"
            subprocess.run([sys.executable, "scripts/validate_submission.py", str(candidate),
                            "--challenge", "9", "--dataset", "testing", "--pairs", str(expected_path),
                            "--report", str(independent_validation)], check=True, stdout=subprocess.DEVNULL)
            with candidate.open() as handle:
                predictions = {tuple(sorted((row["indID1"], row["indID2"]))): row["relation"]
                               for row in csv.DictReader(handle)}
            if set(predictions) != set(truth):
                raise ValueError("Production predictor lost input pairs")
            if baseline_predictions is None:
                baseline_predictions = predictions
            elif predictions != baseline_predictions:
                raise ValueError("Task9 predictions changed under input permutation/reversal")
            summary = json.loads((output / "summary.json").read_text())
            dataset["cases"][case] = {
                "validation": json.loads(independent_validation.read_text()),
                "accuracy": sum(predictions[pair] == label for pair, label in truth.items()) / len(truth),
                "meioses_accuracy": sum(IBD_LABEL_MEIOSES[predictions[pair]] == IBD_LABEL_MEIOSES[label]
                                        for pair, label in truth.items()) / len(truth),
                "prediction_sha256": sha256(candidate), "input_sha256": sha256(observed),
                "changed_predictions": summary["changed_predictions"],
                "release_id": summary["release_id"],
            }
        dataset["permutation_invariant"] = True
        report["datasets"].append(dataset)
    report["passed"] = True
    (args.out_dir / "operational.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
