"""Apply a frozen Task 9 model under the current accepted-label contract."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import joblib
import numpy as np

from ghist_relatedness.ibd import matrix, pair_features, read_segment_rows, read_sex_map
from ghist_relatedness.labels import IBD_LABEL_MEIOSES
from ghist_submission.validate import validate_submission


def accepted_probabilities(classes, probabilities):
    classes = np.asarray(classes)
    keep = np.asarray([str(label) in IBD_LABEL_MEIOSES for label in classes])
    values = np.asarray(probabilities, dtype=float)[:, keep]
    total = values.sum(axis=1, keepdims=True)
    if (
        not values.shape[1]
        or np.any(total <= 0)
        or not np.all(np.isfinite(values))
        or np.any(values < 0)
    ):
        raise ValueError(
            "Model assigned invalid or zero probability to accepted Task 9 labels"
        )
    return classes[keep], values / total


def file_hash(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--observed", type=Path, required=True)
    parser.add_argument("--sex-map", type=Path, required=True)
    parser.add_argument("--dataset", choices=("testing", "final"), required=True)
    parser.add_argument("--evaluation", type=Path)
    parser.add_argument("--genome-cm", type=float, default=3350.933402)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    payload = joblib.load(args.model)
    model, features = payload["model"], payload["features"]
    sex_map = read_sex_map(args.sex_map)
    rows = pair_features(read_segment_rows(args.observed), args.genome_cm, sex_map)
    x, _ = matrix(rows, features)
    raw = model.predict_proba(x)
    labels, probabilities = accepted_probabilities(model.classes_, raw)
    predictions = labels[np.argmax(probabilities, axis=1)]
    original = model.classes_[np.argmax(raw, axis=1)]
    args.out_dir.mkdir(parents=True, exist_ok=True)
    output = args.out_dir / f"task9_{args.dataset}_candidate.csv"
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(("indID1", "indID2", "meioses_count", "relation"))
        for row, label in zip(rows, predictions, strict=True):
            writer.writerow(
                (row["indID1"], row["indID2"], IBD_LABEL_MEIOSES[str(label)], label)
            )
    report = {
        "dataset": args.dataset,
        "policy": "Renormalize the frozen model probabilities over currently accepted labels; no retraining or answer-based selection",
        "changed_predictions": int(np.sum(original != predictions)),
        "validation": validate_submission(
            output, 9, args.dataset, {(r["indID1"], r["indID2"]) for r in rows}
        ),
        "source_sha256": {
            str(p): file_hash(p) for p in (args.model, args.observed, args.sex_map)
        },
        "output_sha256": file_hash(output),
        "submission_count": 0,
    }
    if args.evaluation:
        from scripts.fit_task9_ibd_classifier import load_simulations

        evaluation = load_simulations([args.evaluation], args.genome_cm, sex_map)
        evaluation_x, _ = matrix(evaluation, features)
        evaluation_labels, probability = accepted_probabilities(
            model.classes_, model.predict_proba(evaluation_x)
        )
        predicted = evaluation_labels[np.argmax(probability, axis=1)]
        truth = np.asarray([row["label"] for row in evaluation])
        report["evaluation_on_current_label_set"] = {
            "pairs": len(truth),
            "accuracy": float(np.mean(predicted == truth)),
            "meioses_accuracy": float(
                np.mean(
                    [
                        IBD_LABEL_MEIOSES[str(a)] == IBD_LABEL_MEIOSES[str(b)]
                        for a, b in zip(predicted, truth, strict=True)
                    ]
                )
            ),
            "source": str(args.evaluation),
            "sha256": file_hash(args.evaluation),
        }
    (args.out_dir / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
