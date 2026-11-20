from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ghist_relatedness.validation import validate_simulation_splits
from typing import Any

import joblib
import numpy as np
from sklearn.ensemble import ExtraTreesClassifier, HistGradientBoostingClassifier
from sklearn.metrics import accuracy_score, balanced_accuracy_score
from sklearn.utils.class_weight import compute_sample_weight

try:
    from scripts.fit_human_family_relationships import (
        LABELS,
        MEIOSES,
        feature_names,
        matrix,
        read_observed,
        summarize_raw_lengths,
    )
    from scripts.validate_human_background_correction import (
        read_pedsim_with_lengths,
    )
except ModuleNotFoundError:
    from fit_human_family_relationships import (
        LABELS,
        MEIOSES,
        feature_names,
        matrix,
        read_observed,
        summarize_raw_lengths,
    )
    from validate_human_background_correction import read_pedsim_with_lengths


@dataclass(frozen=True)
class BackgroundRegime:
    mean_segments: float
    excess_scale_cm: float
    pair_heterogeneity: float
    long_tail_probability: float = 0.0
    long_tail_scale: float = 2.0
    minimum_cm: float = 20.0


SELECTION_REGIMES = (
    BackgroundRegime(0.0, 10.0, 0.0),
    BackgroundRegime(1.25, 6.0, 0.10),
    BackgroundRegime(2.50, 14.0, 0.20, 0.05, 2.5),
    BackgroundRegime(3.50, 9.0, 0.30, 0.12, 2.0),
    BackgroundRegime(4.75, 18.0, 0.38, 0.08, 3.0),
    BackgroundRegime(6.00, 24.0, 0.45, 0.15, 2.5),
)

DEVELOPMENT_REGIMES = (
    BackgroundRegime(0.0, 8.0, 0.0),
    BackgroundRegime(1.75, 11.0, 0.15, 0.03, 2.0),
    BackgroundRegime(3.00, 7.0, 0.25, 0.10, 3.0),
    BackgroundRegime(4.00, 16.0, 0.34, 0.06, 2.5),
    BackgroundRegime(5.25, 21.0, 0.42, 0.18, 2.0),
    BackgroundRegime(6.50, 13.0, 0.50, 0.12, 3.5),
)


def random_regime(rng: np.random.Generator) -> BackgroundRegime:
    return BackgroundRegime(
        mean_segments=float(rng.uniform(0.0, 7.0)),
        excess_scale_cm=float(rng.uniform(5.0, 26.0)),
        pair_heterogeneity=float(rng.uniform(0.0, 0.55)),
        long_tail_probability=float(rng.uniform(0.0, 0.20)),
        long_tail_scale=float(rng.uniform(1.5, 3.5)),
    )


def augment_records(
    records: list[dict[str, object]],
    source_lengths: dict[tuple[str, str], list[float]],
    regime: BackgroundRegime,
    rng: np.random.Generator,
    genome_cm: float,
) -> list[dict[str, object]]:
    output = []
    for record in records:
        pair = tuple(sorted((str(record["indID1"]), str(record["indID2"]))))
        multiplier = (
            float(rng.lognormal(-0.5 * regime.pair_heterogeneity**2, regime.pair_heterogeneity))
            if regime.pair_heterogeneity > 0
            else 1.0
        )
        count = int(rng.poisson(regime.mean_segments * multiplier))
        background = regime.minimum_cm + rng.exponential(
            regime.excess_scale_cm, size=count
        )
        if count and regime.long_tail_probability > 0:
            tail = rng.random(count) < regime.long_tail_probability
            background[tail] *= regime.long_tail_scale
        row: dict[str, object] = {
            "indID1": record["indID1"],
            "indID2": record["indID2"],
            "label": record["label"],
        }
        lengths = source_lengths[pair] + background.tolist()
        row.update(summarize_raw_lengths([("0", value) for value in lengths], genome_cm))
        output.append(row)
    return output


def read_replicate(path: Path, genome_cm: float):
    records, lengths = read_pedsim_with_lengths(path, genome_cm)
    return records, lengths


def build_training_set(
    paths: list[Path],
    augmentations: int,
    seed: int,
    genome_cm: float,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    rng = np.random.default_rng(seed)
    output: list[dict[str, object]] = []
    regimes = []
    for path in paths:
        records, lengths = read_replicate(path, genome_cm)
        output.extend(records)
        for _ in range(augmentations):
            regime = random_regime(rng)
            regimes.append(asdict(regime))
            output.extend(augment_records(records, lengths, regime, rng, genome_cm))
    return output, regimes


def regime_datasets(
    path: Path,
    regimes: tuple[BackgroundRegime, ...],
    seed: int,
    genome_cm: float,
) -> list[tuple[BackgroundRegime, list[dict[str, object]]]]:
    records, lengths = read_replicate(path, genome_cm)
    return [
        (
            regime,
            augment_records(
                records,
                lengths,
                regime,
                np.random.default_rng(seed + index),
                genome_cm,
            ),
        )
        for index, regime in enumerate(regimes)
    ]


def fit_direct_models(
    x: np.ndarray, y: np.ndarray, seed: int
) -> list[tuple[str, dict[str, Any]]]:
    models: list[tuple[str, Any]] = [
        (
            "extra_trees_leaf_2",
            ExtraTreesClassifier(
                n_estimators=500,
                min_samples_leaf=2,
                max_features=0.8,
                class_weight="balanced",
                n_jobs=-1,
                random_state=seed,
            ),
        ),
        (
            "extra_trees_leaf_6",
            ExtraTreesClassifier(
                n_estimators=500,
                min_samples_leaf=6,
                max_features=0.8,
                class_weight="balanced",
                n_jobs=-1,
                random_state=seed + 1,
            ),
        ),
        (
            "hist_gradient",
            HistGradientBoostingClassifier(
                learning_rate=0.08,
                max_iter=240,
                max_leaf_nodes=31,
                min_samples_leaf=30,
                l2_regularization=1.0,
                early_stopping=False,
                random_state=seed + 2,
            ),
        ),
    ]
    output = []
    weights = compute_sample_weight("balanced", y)
    for name, model in models:
        if name == "hist_gradient":
            model.fit(x, y, sample_weight=weights)
        else:
            model.fit(x, y)
        output.append((name, {"kind": "direct", "model": model}))
    return output


def fit_hierarchical_model(
    x: np.ndarray, y: np.ndarray, seed: int
) -> tuple[str, dict[str, Any]]:
    degree_y = np.asarray([MEIOSES[str(label)] for label in y])
    degree_model = ExtraTreesClassifier(
        n_estimators=500,
        min_samples_leaf=4,
        max_features=0.8,
        class_weight="balanced",
        n_jobs=-1,
        random_state=seed,
    ).fit(x, degree_y)
    type_models = {}
    for degree in sorted(set(degree_y)):
        selected = degree_y == degree
        labels = np.unique(y[selected])
        if len(labels) == 1:
            type_models[int(degree)] = {"constant": str(labels[0])}
            continue
        model = ExtraTreesClassifier(
            n_estimators=400,
            min_samples_leaf=4,
            max_features=1.0,
            class_weight="balanced",
            n_jobs=-1,
            random_state=seed + int(degree),
        ).fit(x[selected], y[selected])
        type_models[int(degree)] = {"model": model}
    return (
        "hierarchical_extra_trees",
        {"kind": "hierarchical", "degree_model": degree_model, "type_models": type_models},
    )


def predict_probabilities(payload: dict[str, Any], x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    if payload["kind"] == "direct":
        model = payload["model"]
        return np.asarray(model.classes_), model.predict_proba(x)
    degree_model = payload["degree_model"]
    degree_probabilities = degree_model.predict_proba(x)
    classes = np.asarray(LABELS)
    probabilities = np.zeros((x.shape[0], len(classes)), dtype=float)
    class_index = {str(label): index for index, label in enumerate(classes)}
    for degree_index, degree in enumerate(degree_model.classes_):
        type_payload = payload["type_models"][int(degree)]
        if "constant" in type_payload:
            probabilities[:, class_index[type_payload["constant"]]] = degree_probabilities[:, degree_index]
            continue
        type_model = type_payload["model"]
        conditional = type_model.predict_proba(x)
        for type_index, label in enumerate(type_model.classes_):
            probabilities[:, class_index[str(label)]] = (
                degree_probabilities[:, degree_index] * conditional[:, type_index]
            )
    return classes, probabilities


def metrics(payload: dict[str, Any], records: list[dict[str, object]], names: list[str]):
    truth = np.asarray([str(record["label"]) for record in records])
    classes, probabilities = predict_probabilities(payload, matrix(records, names))
    predicted = classes[np.argmax(probabilities, axis=1)]
    truth_degree = np.asarray([MEIOSES[label] for label in truth])
    predicted_degree = np.asarray([MEIOSES[str(label)] for label in predicted])
    return {
        "exact_accuracy": float(accuracy_score(truth, predicted)),
        "balanced_accuracy": float(balanced_accuracy_score(truth, predicted)),
        "meioses_accuracy": float(accuracy_score(truth_degree, predicted_degree)),
    }


def evaluate_regimes(
    payload: dict[str, Any],
    datasets: list[tuple[BackgroundRegime, list[dict[str, object]]]],
    names: list[str],
) -> dict[str, object]:
    runs = []
    for regime, records in datasets:
        row = asdict(regime)
        row.update(metrics(payload, records, names))
        runs.append(row)
    return {
        "mean_exact_accuracy": float(np.mean([row["exact_accuracy"] for row in runs])),
        "worst_exact_accuracy": float(np.min([row["exact_accuracy"] for row in runs])),
        "mean_balanced_accuracy": float(np.mean([row["balanced_accuracy"] for row in runs])),
        "mean_meioses_accuracy": float(np.mean([row["meioses_accuracy"] for row in runs])),
        "runs": runs,
    }


def observed_predictions(
    payload: dict[str, Any], observed: list[dict[str, object]], names: list[str]
) -> list[dict[str, object]]:
    classes, probabilities = predict_probabilities(payload, matrix(observed, names))
    predictions = classes[np.argmax(probabilities, axis=1)]
    return [
        {
            "indID1": record["indID1"],
            "indID2": record["indID2"],
            "relation": str(predictions[index]),
            "meioses_count": MEIOSES[str(predictions[index])],
            "confidence": float(np.max(probabilities[index])),
            "probabilities": {
                str(label): float(probability)
                for label, probability in zip(classes, probabilities[index], strict=True)
            },
        }
        for index, record in enumerate(observed)
    ]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--calibration", type=Path, nargs="+", required=True)
    parser.add_argument("--selection", type=Path, required=True)
    parser.add_argument("--development-evaluation", type=Path, required=True)
    parser.add_argument("--observed", type=Path, required=True)
    parser.add_argument("--minimum-cm", type=int, default=20)
    parser.add_argument("--augmentations", type=int, default=4)
    parser.add_argument("--genome-cm", type=float, default=3346.298892538)
    parser.add_argument("--seed", type=int, default=120826)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    validate_simulation_splits({
        'calibration': args.calibration,
        'selection': [args.selection],
        'development_evaluation': [args.development_evaluation],
    })

    names = feature_names(args.minimum_cm)
    train, training_regimes = build_training_set(
        args.calibration,
        args.augmentations,
        args.seed,
        args.genome_cm,
    )
    train_x = matrix(train, names)
    train_y = np.asarray([str(record["label"]) for record in train])
    candidates = fit_direct_models(train_x, train_y, args.seed)
    candidates.append(fit_hierarchical_model(train_x, train_y, args.seed + 100))

    selection = regime_datasets(
        args.selection, SELECTION_REGIMES, args.seed + 1000, args.genome_cm
    )
    results = []
    fitted = {}
    for name, payload in candidates:
        result = evaluate_regimes(payload, selection, names)
        result["name"] = name
        results.append(result)
        fitted[name] = payload
    results.sort(
        key=lambda row: (
            row["mean_exact_accuracy"],
            row["worst_exact_accuracy"],
            row["mean_meioses_accuracy"],
        ),
        reverse=True,
    )
    selected_name = str(results[0]["name"])
    selected = fitted[selected_name]
    development = regime_datasets(
        args.development_evaluation,
        DEVELOPMENT_REGIMES,
        args.seed + 2000,
        args.genome_cm,
    )
    development_result = evaluate_regimes(selected, development, names)
    observed = read_observed(args.observed)
    predictions = observed_predictions(selected, observed, names)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    joblib.dump(
        {"payload": selected, "feature_names": names}, args.out_dir / "model.joblib"
    )
    (args.out_dir / "predictions.json").write_text(
        json.dumps(predictions, indent=2) + "\n"
    )
    with (args.out_dir / "candidate.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("indID1", "indID2", "relation"))
        writer.writeheader()
        writer.writerows(
            {key: row[key] for key in writer.fieldnames} for row in predictions
        )
    summary = {
        "state": "DEVELOPMENT",
        "design": "background-augmented segment-shape inference with a hierarchical degree/type candidate",
        "selected_model": selected_name,
        "minimum_cm": args.minimum_cm,
        "training_pairs": len(train),
        "training_regimes": training_regimes,
        "selection": results,
        "development_evaluation": development_result,
        "observed_pairs": len(predictions),
        "observed_label_counts": {
            label: sum(row["relation"] == label for row in predictions)
            for label in LABELS
        },
        "final_holdout_evaluated": False,
        "submission_count": 0,
    }
    (args.out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({key: value for key, value in summary.items() if key != "selection"}, indent=2))


if __name__ == "__main__":
    main()
