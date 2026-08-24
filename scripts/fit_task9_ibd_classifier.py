from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import joblib
import numpy as np
from sklearn.ensemble import ExtraTreesClassifier, HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.metrics import accuracy_score, balanced_accuracy_score, confusion_matrix

from ghist_relatedness.ibd import (
    SexMap,
    matrix,
    pair_features,
    read_segment_rows,
    read_sex_map,
    relationship_from_sample,
)


MEIOSES = {
    "siblings": 2,
    "half_sibling": 2,
    "grandparent_grandchild": 2,
    "avuncular": 3,
    "double_cousin": 4,
    "full_1st_cousin": 4,
    "half_1st_cousin": 4,
    "full_2nd_cousin": 6,
    "half_2nd_cousin": 6,
    "full_3rd_cousin": 8,
    "full_4th_cousin": 10,
    "full_5th_cousin": 12,
}


def load_simulations(
    paths: list[Path], genome_cm: float, sex_map: SexMap
) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for path in paths:
        for record in pair_features(read_segment_rows(path), genome_cm, sex_map):
            label1 = relationship_from_sample(str(record["indID1"]))
            label2 = relationship_from_sample(str(record["indID2"]))
            if label1 != label2:
                raise ValueError(f"Mismatched simulated pair in {path}")
            record["label"] = label1
            records.append(record)
    return records


def score(model, records: list[dict[str, object]], feature_names: list[str]) -> dict[str, object]:
    features, _ = matrix(records, feature_names)
    truth = np.asarray([record["label"] for record in records])
    predicted = model.predict(features)
    labels = sorted(MEIOSES)
    return {
        "pairs": len(records),
        "accuracy": float(accuracy_score(truth, predicted)),
        "balanced_accuracy": float(balanced_accuracy_score(truth, predicted)),
        "labels": labels,
        "confusion_matrix": confusion_matrix(truth, predicted, labels=labels).tolist(),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--calibration", type=Path, nargs="+", required=True)
    parser.add_argument("--validation", type=Path, required=True)
    parser.add_argument("--evaluation", type=Path, required=True)
    parser.add_argument("--observed", type=Path, required=True)
    parser.add_argument("--sex-map", type=Path, required=True)
    parser.add_argument("--genome-cm", type=float, default=3350.933402)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()

    sex_map = read_sex_map(args.sex_map)
    calibration = load_simulations(args.calibration, args.genome_cm, sex_map)
    validation = load_simulations([args.validation], args.genome_cm, sex_map)
    evaluation = load_simulations([args.evaluation], args.genome_cm, sex_map)
    train_x, feature_names = matrix(calibration)
    train_y = np.asarray([record["label"] for record in calibration])
    candidates = {
        "extra_trees": ExtraTreesClassifier(
            n_estimators=600,
            min_samples_leaf=2,
            max_features=0.8,
            class_weight="balanced",
            n_jobs=-1,
            random_state=90926,
        ),
        "random_forest": RandomForestClassifier(
            n_estimators=600,
            min_samples_leaf=2,
            max_features=0.8,
            class_weight="balanced_subsample",
            n_jobs=-1,
            random_state=90926,
        ),
        "hist_gradient_boosting": HistGradientBoostingClassifier(
            max_iter=300,
            learning_rate=0.08,
            max_leaf_nodes=31,
            l2_regularization=1.0,
            random_state=90926,
        ),
    }
    validation_scores = {}
    for name, model in candidates.items():
        model.fit(train_x, train_y)
        validation_scores[name] = score(model, validation, feature_names)
    selected_name = max(
        candidates,
        key=lambda name: (
            validation_scores[name]["balanced_accuracy"],
            validation_scores[name]["accuracy"],
            name == "extra_trees",
        ),
    )
    selected = candidates[selected_name]
    evaluation_score = score(selected, evaluation, feature_names)

    observed_records = pair_features(read_segment_rows(args.observed), args.genome_cm, sex_map)
    observed_x, _ = matrix(observed_records, feature_names)
    probabilities = selected.predict_proba(observed_x)
    predictions = selected.classes_[np.argmax(probabilities, axis=1)]
    confidence = np.max(probabilities, axis=1)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    prediction_path = args.out_dir / "task9_testing_predictions.csv"
    with prediction_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(("indID1", "indID2", "relationship", "meioses", "confidence"))
        for record, label, probability in zip(observed_records, predictions, confidence, strict=True):
            writer.writerow((record["indID1"], record["indID2"], label, MEIOSES[str(label)], f"{probability:.8f}"))

    joblib.dump({"model": selected, "features": feature_names}, args.out_dir / "task9_model.joblib")
    summary = {
        "state": "COMPLETED",
        "protocol": {
            "calibration_replicates": [str(path) for path in args.calibration],
            "validation_replicate": str(args.validation),
            "evaluation_replicate": str(args.evaluation),
            "selection_rule": "maximum validation balanced accuracy; evaluation untouched until selection",
        },
        "feature_names": feature_names,
        "validation_candidates": validation_scores,
        "selected_model": selected_name,
        "locked_evaluation": evaluation_score,
        "observed_pairs": len(observed_records),
        "observed_label_counts": {
            label: int(np.sum(predictions == label)) for label in selected.classes_
        },
        "observed_confidence": {
            "median": float(np.median(confidence)),
            "q10": float(np.quantile(confidence, 0.10)),
            "minimum": float(np.min(confidence)),
        },
        "prediction_path": str(prediction_path),
        "submission_count": 0,
    }
    (args.out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
