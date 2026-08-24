from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

import joblib
import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.utils.class_weight import compute_sample_weight

try:
    from scripts.fit_human_background_corrected import align_probabilities
    from scripts.fit_human_family_relationships import (
        LABELS,
        MEIOSES,
        feature_names,
        matrix,
        read_observed,
        read_pedsim,
    )
except ModuleNotFoundError:
    from fit_human_background_corrected import align_probabilities
    from fit_human_family_relationships import (
        LABELS,
        MEIOSES,
        feature_names,
        matrix,
        read_observed,
        read_pedsim,
    )


def spectrum_feature_names(minimum_cm: int) -> list[str]:
    return [
        name
        for name in feature_names(minimum_cm)
        if name != "max_cm" and not name.startswith("ordered_length_")
    ]


def model_specs():
    for leaf_nodes, minimum_leaf, regularization in (
        (15, 20, 1.0),
        (31, 20, 1.0),
        (31, 40, 3.0),
        (63, 30, 5.0),
    ):
        yield {
            "max_leaf_nodes": leaf_nodes,
            "min_samples_leaf": minimum_leaf,
            "l2_regularization": regularization,
        }


def fit_models(
    records: list[dict[str, object]], minimum_cm: int, seed: int
):
    labels = np.asarray([str(record["label"]) for record in records])
    weights = compute_sample_weight("balanced", labels)
    feature_sets = {
        "full_shape": feature_names(minimum_cm),
        "segment_spectrum": spectrum_feature_names(minimum_cm),
    }
    fitted = []
    for feature_set, names in feature_sets.items():
        values = matrix(records, names)
        for index, spec in enumerate(model_specs()):
            model = HistGradientBoostingClassifier(
                learning_rate=0.07,
                max_iter=300,
                early_stopping=True,
                random_state=seed + 10 * len(fitted) + index,
                **spec,
            ).fit(values, labels, sample_weight=weights)
            fitted.append(
                {
                    "name": f"{feature_set}_leaves_{spec['max_leaf_nodes']}_min_{spec['min_samples_leaf']}_l2_{spec['l2_regularization']}",
                    "feature_set": feature_set,
                    "feature_names": names,
                    "model": model,
                    "parameters": spec,
                }
            )
    return fitted


def metrics(classes, probabilities, records):
    truth = np.asarray([str(record["label"]) for record in records])
    predicted = classes[np.argmax(probabilities, axis=1)]
    return {
        "exact_accuracy": float(np.mean(predicted == truth)),
        "meioses_accuracy": float(
            np.mean(
                [
                    MEIOSES[str(left)] == MEIOSES[str(right)]
                    for left, right in zip(predicted, truth, strict=True)
                ]
            )
        ),
    }


def blended_probabilities(candidate, baseline, records, weight):
    classes = np.asarray(LABELS)
    candidate_probabilities = candidate["model"].predict_proba(
        matrix(records, candidate["feature_names"])
    )
    baseline_probabilities = baseline["model"].predict_proba(
        matrix(records, baseline["feature_names"])
    )
    return classes, weight * align_probabilities(
        np.asarray(candidate["model"].classes_), candidate_probabilities, classes
    ) + (1.0 - weight) * align_probabilities(
        np.asarray(baseline["model"].classes_), baseline_probabilities, classes
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline-model", type=Path, required=True)
    parser.add_argument("--calibration", type=Path, nargs="+", required=True)
    parser.add_argument("--selection", type=Path, required=True)
    parser.add_argument("--development-evaluation", type=Path, required=True)
    parser.add_argument("--observed", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--minimum-cm", type=int, default=5)
    parser.add_argument("--genome-cm", type=float, default=3346.298892538)
    parser.add_argument("--seed", type=int, default=102426)
    args = parser.parse_args()

    baseline = joblib.load(args.baseline_model)
    training = [
        record
        for path in args.calibration
        for record in read_pedsim(path, args.genome_cm)
    ]
    selection = read_pedsim(args.selection, args.genome_cm)
    development = read_pedsim(args.development_evaluation, args.genome_cm)
    fitted = fit_models(training, args.minimum_cm, args.seed)
    results = []
    for candidate in fitted:
        for weight in (0.25, 0.50, 0.75, 1.0):
            classes, probabilities = blended_probabilities(
                candidate, baseline, selection, weight
            )
            result = metrics(classes, probabilities, selection)
            result.update(
                {
                    "name": candidate["name"],
                    "candidate_weight": weight,
                    "feature_set": candidate["feature_set"],
                    "parameters": candidate["parameters"],
                }
            )
            results.append(result)
    results.sort(
        key=lambda row: (row["exact_accuracy"], row["meioses_accuracy"]),
        reverse=True,
    )
    winner = results[0]
    selected = next(item for item in fitted if item["name"] == winner["name"])
    weight = float(winner["candidate_weight"])
    classes, probabilities = blended_probabilities(
        selected, baseline, development, weight
    )
    development_result = metrics(classes, probabilities, development)
    baseline_classes = np.asarray(baseline["model"].classes_)
    baseline_development = metrics(
        baseline_classes,
        baseline["model"].predict_proba(
            matrix(development, baseline["feature_names"])
        ),
        development,
    )
    observed = read_observed(args.observed)
    classes, probabilities = blended_probabilities(selected, baseline, observed, weight)
    predicted = classes[np.argmax(probabilities, axis=1)]
    rows = [
        {
            "indID1": record["indID1"],
            "indID2": record["indID2"],
            "relation": str(predicted[index]),
            "confidence": float(np.max(probabilities[index])),
            "probabilities": {
                str(label): float(value)
                for label, value in zip(classes, probabilities[index], strict=True)
            },
        }
        for index, record in enumerate(observed)
    ]
    args.out_dir.mkdir(parents=True, exist_ok=True)
    joblib.dump(
        {
            "candidate_model": selected["model"],
            "candidate_feature_names": selected["feature_names"],
            "candidate_weight": weight,
        },
        args.out_dir / "model.joblib",
    )
    (args.out_dir / "predictions.json").write_text(json.dumps(rows, indent=2) + "\n")
    with (args.out_dir / "candidate.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("indID1", "indID2", "relation"))
        writer.writeheader()
        writer.writerows({key: row[key] for key in writer.fieldnames} for row in rows)
    summary = {
        "state": "DEVELOPMENT",
        "design": "histogram-gradient comparison of full ordered shape and compact segment-spectrum features",
        "training_pairs": len(training),
        "selection_pairs": len(selection),
        "development_pairs": len(development),
        "selection": results,
        "selected_model": winner,
        "development_evaluation": development_result,
        "development_baseline": baseline_development,
        "observed_pairs": len(rows),
        "observed_label_counts": dict(Counter(row["relation"] for row in rows)),
        "final_holdout_evaluated": False,
        "submission_count": 0,
    }
    (args.out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(
        json.dumps(
            {
                "selected_model": winner,
                "development_evaluation": development_result,
                "development_baseline": baseline_development,
                "observed_label_counts": summary["observed_label_counts"],
                "submission_count": 0,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
