from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

from ghist_relatedness.validation import validate_simulation_splits

import joblib
import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import balanced_accuracy_score
from sklearn.utils.class_weight import compute_sample_weight

from ghist_relatedness.ibd import matrix, pair_features, read_segment_rows, read_sex_map
from ghist_relatedness.labels import IBD_LABEL_MEIOSES
try:
    from scripts.fit_task9_ibd_classifier import MEIOSES, load_simulations
except ModuleNotFoundError:
    from fit_task9_ibd_classifier import MEIOSES, load_simulations


LABELS = tuple(sorted(IBD_LABEL_MEIOSES))


def classifier(seed: int, leaves: int, minimum_leaf: int):
    return HistGradientBoostingClassifier(
        learning_rate=0.07,
        max_iter=300,
        max_leaf_nodes=leaves,
        min_samples_leaf=minimum_leaf,
        l2_regularization=3.0,
        early_stopping=False,
        random_state=seed,
    )


def aligned(classes, probabilities):
    output = np.zeros((len(probabilities), len(LABELS)), dtype=float)
    indices = {label: index for index, label in enumerate(LABELS)}
    for source, label in enumerate(classes):
        if str(label) in indices:
            output[:, indices[str(label)]] = probabilities[:, source]
    totals = output.sum(axis=1, keepdims=True)
    if np.any(totals <= 0):
        raise ValueError("Model assigned no probability to accepted Task 9 labels")
    return output / totals


def fit_models(x, y, seed, leaves, minimum_leaf):
    exact = classifier(seed, leaves, minimum_leaf).fit(
        x, y, sample_weight=compute_sample_weight("balanced", y)
    )
    degrees = np.asarray([MEIOSES[str(label)] for label in y])
    degree = classifier(seed + 1, leaves, minimum_leaf).fit(
        x, degrees, sample_weight=compute_sample_weight("balanced", degrees)
    )
    conditional = {}
    for meioses in sorted(set(degrees)):
        selected = degrees == meioses
        if len(np.unique(y[selected])) > 1:
            conditional[int(meioses)] = classifier(
                seed + 10 + int(meioses), leaves, minimum_leaf
            ).fit(
                x[selected],
                y[selected],
                sample_weight=compute_sample_weight("balanced", y[selected]),
            )
    return {"exact": exact, "degree": degree, "conditional": conditional}


def factorized_probabilities(bundle, x):
    degree_probabilities = bundle["degree"].predict_proba(x)
    output = np.zeros((len(x), len(LABELS)), dtype=float)
    indices = {label: index for index, label in enumerate(LABELS)}
    for column, meioses in enumerate(bundle["degree"].classes_):
        labels = [str(label) for label in bundle["exact"].classes_ if MEIOSES[str(label)] == int(meioses)]
        if len(labels) == 1:
            output[:, indices[labels[0]]] = degree_probabilities[:, column]
        else:
            model = bundle["conditional"][int(meioses)]
            conditional = model.predict_proba(x)
            for source, label in enumerate(model.classes_):
                output[:, indices[str(label)]] = degree_probabilities[:, column] * conditional[:, source]
    return output / np.maximum(output.sum(axis=1, keepdims=True), 1e-12)


def metrics(probabilities, truth):
    predicted = np.asarray(LABELS)[np.argmax(probabilities, axis=1)]
    return {
        "accuracy": float(np.mean(predicted == truth)),
        "balanced_accuracy": float(balanced_accuracy_score(truth, predicted)),
        "meioses_accuracy": float(
            np.mean([MEIOSES[str(a)] == MEIOSES[str(b)] for a, b in zip(predicted, truth, strict=True)])
        ),
    }


def combined_probabilities(bundle, baseline, x, winner):
    direct = aligned(bundle["exact"].classes_, bundle["exact"].predict_proba(x))
    factorized = factorized_probabilities(bundle, x)
    structured = winner["factor_weight"] * factorized + (1.0 - winner["factor_weight"]) * direct
    baseline_probabilities = aligned(baseline.classes_, baseline.predict_proba(x))
    return winner["candidate_weight"] * structured + (1.0 - winner["candidate_weight"]) * baseline_probabilities


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline-model", type=Path, required=True)
    parser.add_argument("--calibration", type=Path, nargs="+", required=True)
    parser.add_argument("--selection", type=Path, required=True)
    parser.add_argument("--development-evaluation", type=Path, required=True)
    parser.add_argument("--observed", type=Path, required=True)
    parser.add_argument("--sex-map", type=Path, required=True)
    parser.add_argument("--genome-cm", type=float, default=3350.933402)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=92026)
    args = parser.parse_args()
    validate_simulation_splits({
        'calibration': args.calibration,
        'selection': [args.selection],
        'development_evaluation': [args.development_evaluation],
    })

    sex_map = read_sex_map(args.sex_map)
    calibration = load_simulations(args.calibration, args.genome_cm, sex_map)
    selection = load_simulations([args.selection], args.genome_cm, sex_map)
    development = load_simulations([args.development_evaluation], args.genome_cm, sex_map)
    train_x, features = matrix(calibration)
    selection_x, _ = matrix(selection, features)
    development_x, _ = matrix(development, features)
    train_y = np.asarray([row["label"] for row in calibration])
    selection_y = np.asarray([row["label"] for row in selection])
    development_y = np.asarray([row["label"] for row in development])
    baseline_payload = joblib.load(args.baseline_model)
    baseline = baseline_payload["model"]
    baseline_selection = aligned(baseline.classes_, baseline.predict_proba(selection_x))
    baseline_selection_metrics = metrics(baseline_selection, selection_y)

    screens = []
    fitted = []
    for leaves, minimum_leaf in ((15, 30), (31, 30), (31, 60)):
        bundle = fit_models(train_x, train_y, args.seed + 100 * len(fitted), leaves, minimum_leaf)
        direct = aligned(bundle["exact"].classes_, bundle["exact"].predict_proba(selection_x))
        factorized = factorized_probabilities(bundle, selection_x)
        for factor_weight in (0.25, 0.50, 0.75):
            structured = factor_weight * factorized + (1.0 - factor_weight) * direct
            for candidate_weight in (0.25, 0.50, 0.75):
                probabilities = candidate_weight * structured + (1.0 - candidate_weight) * baseline_selection
                result = metrics(probabilities, selection_y)
                result.update(
                    {
                        "leaves": leaves,
                        "minimum_leaf": minimum_leaf,
                        "factor_weight": factor_weight,
                        "candidate_weight": candidate_weight,
                    }
                )
                screens.append(result)
        fitted.append((leaves, minimum_leaf, bundle))
    eligible = [
        row
        for row in screens
        if row["meioses_accuracy"] >= baseline_selection_metrics["meioses_accuracy"] - 0.001
        and row["balanced_accuracy"] >= baseline_selection_metrics["balanced_accuracy"]
    ]
    eligible.sort(key=lambda row: (row["accuracy"], row["balanced_accuracy"], row["meioses_accuracy"]), reverse=True)
    if not eligible:
        raise ValueError("No candidate passed the selection non-regression guard")
    winner = eligible[0]
    bundle = next(item[2] for item in fitted if item[:2] == (winner["leaves"], winner["minimum_leaf"]))
    development_probabilities = combined_probabilities(bundle, baseline, development_x, winner)
    baseline_development = aligned(baseline.classes_, baseline.predict_proba(development_x))

    observed_records = pair_features(read_segment_rows(args.observed), args.genome_cm, sex_map)
    observed_x, _ = matrix(observed_records, features)
    observed_probabilities = combined_probabilities(bundle, baseline, observed_x, winner)
    observed_labels = np.asarray(LABELS)[np.argmax(observed_probabilities, axis=1)]
    rows = [
        {
            "indID1": record["indID1"],
            "indID2": record["indID2"],
            "relation": str(observed_labels[index]),
            "probabilities": {label: float(value) for label, value in zip(LABELS, observed_probabilities[index], strict=True)},
        }
        for index, record in enumerate(observed_records)
    ]
    args.out_dir.mkdir(parents=True, exist_ok=True)
    joblib.dump({"bundle": bundle, "winner": winner, "features": features}, args.out_dir / "model.joblib")
    (args.out_dir / "predictions.json").write_text(json.dumps(rows, indent=2) + "\n")
    with (args.out_dir / "task9_candidate.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(("indID1", "indID2", "meioses_count", "relation"))
        for row in rows:
            writer.writerow((row["indID1"], row["indID2"], MEIOSES[row["relation"]], row["relation"]))
    summary = {
        "state": "DEVELOPMENT",
        "design": "degree-factorized histogram-gradient and calibrated-forest blend",
        "selected": winner,
        "selection_baseline": baseline_selection_metrics,
        "selection": screens,
        "development_evaluation": metrics(development_probabilities, development_y),
        "development_baseline": metrics(baseline_development, development_y),
        "observed_pairs": len(rows),
        "observed_label_counts": dict(Counter(observed_labels)),
        "submission_count": 0,
    }
    (args.out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({key: value for key, value in summary.items() if key != "selection"}, indent=2))


if __name__ == "__main__":
    main()
