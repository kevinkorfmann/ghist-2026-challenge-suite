from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

import joblib
import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.utils.class_weight import compute_sample_weight

try:
    from scripts.fit_task11_dog_relationships import (
        LABELS,
        feature_names,
        read_clusters,
        read_observed,
        relationship_label,
        score,
        summarize_segments,
    )
except ModuleNotFoundError:
    from fit_task11_dog_relationships import (
        LABELS,
        feature_names,
        read_clusters,
        read_observed,
        relationship_label,
        score,
        summarize_segments,
    )


MEIOSES = {
    "siblings": 2,
    "half_sibling": 2,
    "avuncular": 3,
    "half_avuncular": 3,
    "avuncular_once_removed": 4,
    "full_1st_cousin": 4,
    "half_1st_cousin": 4,
    "avuncular_twice_removed": 5,
}
STATE_TOP = 8


def rich_feature_names() -> list[str]:
    names = feature_names()
    for state in ("ibd1", "ibd2"):
        names.extend(
            [
                f"{state}_segments",
                f"{state}_chromosomes",
                f"{state}_mean",
                f"{state}_median",
                f"{state}_q75",
                f"{state}_q90",
                f"{state}_max",
            ]
        )
        names.extend(f"{state}_ordered_{rank}" for rank in range(1, STATE_TOP + 1))
    names.extend(
        [
            "chromosome_sharing_mean",
            "chromosome_sharing_std",
            "chromosome_sharing_q75",
            "chromosome_sharing_max",
            "chromosomes_with_ibd2",
        ]
    )
    return names


def rich_summary(
    records: list[tuple[str, float, bool]], exposure_scale: float = 1.0
) -> np.ndarray:
    values = summarize_segments(records, exposure_scale)
    extensions: list[float] = []
    for full in (False, True):
        selected = [record for record in records if record[2] is full]
        lengths = np.asarray([record[1] for record in selected], dtype=float)
        extensions.extend([float(len(selected)), float(len({row[0] for row in selected}))])
        if lengths.size:
            extensions.extend(
                [
                    float(lengths.mean()),
                    float(np.median(lengths)),
                    float(np.quantile(lengths, 0.75)),
                    float(np.quantile(lengths, 0.90)),
                    float(lengths.max()),
                ]
            )
            ordered = np.sort(lengths)[::-1][:STATE_TOP]
        else:
            extensions.extend([0.0] * 5)
            ordered = np.asarray([], dtype=float)
        extensions.extend(np.pad(ordered, (0, STATE_TOP - ordered.size)).tolist())

    by_chromosome: dict[str, float] = defaultdict(float)
    chromosomes_with_ibd2 = set()
    for chromosome, length, full in records:
        by_chromosome[chromosome] += length * (2.0 if full else 1.0) * exposure_scale
        if full:
            chromosomes_with_ibd2.add(chromosome)
    sharing = np.asarray(list(by_chromosome.values()), dtype=float)
    extensions.extend(
        [
            float(sharing.mean()) if sharing.size else 0.0,
            float(sharing.std()) if sharing.size else 0.0,
            float(np.quantile(sharing, 0.75)) if sharing.size else 0.0,
            float(sharing.max()) if sharing.size else 0.0,
            float(len(chromosomes_with_ibd2)),
        ]
    )
    return np.concatenate((values, np.asarray(extensions)))


def read_pedsim_rich(path: Path) -> tuple[np.ndarray, np.ndarray]:
    grouped: dict[tuple[str, str], list[tuple[str, float, bool]]] = defaultdict(list)
    labels: dict[tuple[str, str], str] = {}
    with path.open() as handle:
        for line in handle:
            left, right, chromosome, _, _, state, _, _, length = line.split()
            pair = tuple(sorted((left, right)))
            label = relationship_label(left)
            if relationship_label(right) != label:
                raise ValueError(f"Mixed relationship labels for {pair}")
            labels[pair] = label
            grouped[pair].append((chromosome, float(length), state == "IBD2"))
    pairs = sorted(grouped)
    return (
        np.vstack([rich_summary(grouped[pair]) for pair in pairs]),
        np.asarray([labels[pair] for pair in pairs]),
    )


def aligned(classes: np.ndarray, probabilities: np.ndarray) -> np.ndarray:
    output = np.zeros((len(probabilities), len(LABELS)), dtype=float)
    indices = {label: index for index, label in enumerate(LABELS)}
    for source, label in enumerate(classes):
        output[:, indices[str(label)]] = probabilities[:, source]
    return output


def new_classifier(seed: int, max_leaf_nodes: int, minimum_leaf: int):
    return HistGradientBoostingClassifier(
        learning_rate=0.07,
        max_iter=300,
        max_leaf_nodes=max_leaf_nodes,
        min_samples_leaf=minimum_leaf,
        l2_regularization=3.0,
        early_stopping=True,
        random_state=seed,
    )


def fit_factorized(x: np.ndarray, y: np.ndarray, seed: int, leaves: int, minimum_leaf: int):
    weights = compute_sample_weight("balanced", y)
    exact = new_classifier(seed, leaves, minimum_leaf).fit(x, y, sample_weight=weights)
    degrees = np.asarray([MEIOSES[str(label)] for label in y])
    degree = new_classifier(seed + 1, leaves, minimum_leaf).fit(
        x, degrees, sample_weight=compute_sample_weight("balanced", degrees)
    )
    conditional = {}
    for meioses in sorted(set(degrees)):
        selected = degrees == meioses
        labels = np.unique(y[selected])
        if len(labels) > 1:
            conditional[int(meioses)] = new_classifier(
                seed + 10 + int(meioses), leaves, minimum_leaf
            ).fit(
                x[selected],
                y[selected],
                sample_weight=compute_sample_weight("balanced", y[selected]),
            )
    return {"exact": exact, "degree": degree, "conditional": conditional}


def factorized_probabilities(bundle, x: np.ndarray) -> np.ndarray:
    degree_probabilities = bundle["degree"].predict_proba(x)
    output = np.zeros((len(x), len(LABELS)), dtype=float)
    label_index = {label: index for index, label in enumerate(LABELS)}
    for column, meioses in enumerate(bundle["degree"].classes_):
        labels = [label for label in LABELS if MEIOSES[label] == int(meioses)]
        if len(labels) == 1:
            output[:, label_index[labels[0]]] = degree_probabilities[:, column]
            continue
        model = bundle["conditional"][int(meioses)]
        conditional = model.predict_proba(x)
        for source, label in enumerate(model.classes_):
            output[:, label_index[str(label)]] = degree_probabilities[:, column] * conditional[:, source]
    return output / np.maximum(output.sum(axis=1, keepdims=True), 1e-12)


def metrics(probabilities: np.ndarray, truth: np.ndarray) -> dict[str, float]:
    predicted = np.asarray(LABELS)[np.argmax(probabilities, axis=1)]
    return {
        "accuracy": float(np.mean(predicted == truth)),
        "meioses_accuracy": float(
            np.mean([MEIOSES[str(a)] == MEIOSES[str(b)] for a, b in zip(predicted, truth, strict=True)])
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline-model", type=Path, required=True)
    parser.add_argument("--calibration", type=Path, nargs="+", required=True)
    parser.add_argument("--selection", type=Path, required=True)
    parser.add_argument("--development-evaluation", type=Path, required=True)
    parser.add_argument("--observed-ibd", type=Path, required=True)
    parser.add_argument("--observed-hbd", type=Path, required=True)
    parser.add_argument("--hbd-threshold-cm", type=float, default=5.0)
    parser.add_argument("--map", type=Path, required=True)
    parser.add_argument("--clusters", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=112026)
    args = parser.parse_args()

    training = [read_pedsim_rich(path) for path in args.calibration]
    train_x = np.vstack([item[0] for item in training])
    train_y = np.concatenate([item[1] for item in training])
    selection_x, selection_y = read_pedsim_rich(args.selection)
    development_x, development_y = read_pedsim_rich(args.development_evaluation)
    baseline = joblib.load(args.baseline_model)

    screens = []
    fitted = []
    for leaves, minimum_leaf in ((15, 30), (31, 30), (31, 60)):
        bundle = fit_factorized(train_x, train_y, args.seed + len(fitted) * 100, leaves, minimum_leaf)
        factorized = factorized_probabilities(bundle, selection_x)
        direct = aligned(np.asarray(bundle["exact"].classes_), bundle["exact"].predict_proba(selection_x))
        baseline_probabilities = aligned(np.asarray(baseline.classes_), baseline.predict_proba(selection_x[:, : len(feature_names())]))
        for factor_weight in (0.5, 0.75, 1.0):
            structured = factor_weight * factorized + (1.0 - factor_weight) * direct
            for candidate_weight in (0.5, 0.75, 1.0):
                probabilities = candidate_weight * structured + (1.0 - candidate_weight) * baseline_probabilities
                result = metrics(probabilities, selection_y)
                result.update({"leaves": leaves, "minimum_leaf": minimum_leaf, "factor_weight": factor_weight, "candidate_weight": candidate_weight})
                screens.append(result)
        fitted.append((leaves, minimum_leaf, bundle))
    screens.sort(key=lambda row: (row["accuracy"], row["meioses_accuracy"]), reverse=True)
    winner = screens[0]
    bundle = next(item[2] for item in fitted if item[:2] == (winner["leaves"], winner["minimum_leaf"]))

    factorized = factorized_probabilities(bundle, development_x)
    direct = aligned(np.asarray(bundle["exact"].classes_), bundle["exact"].predict_proba(development_x))
    baseline_development = aligned(np.asarray(baseline.classes_), baseline.predict_proba(development_x[:, : len(feature_names())]))
    structured = winner["factor_weight"] * factorized + (1.0 - winner["factor_weight"]) * direct
    development = winner["candidate_weight"] * structured + (1.0 - winner["candidate_weight"]) * baseline_development

    observed_x, observed_pairs, family_indices = read_observed(
        args.observed_ibd, args.observed_hbd, args.map, args.clusters, args.hbd_threshold_cm, rich_summary
    )
    factorized = factorized_probabilities(bundle, observed_x)
    direct = aligned(np.asarray(bundle["exact"].classes_), bundle["exact"].predict_proba(observed_x))
    baseline_observed = aligned(np.asarray(baseline.classes_), baseline.predict_proba(observed_x[:, : len(feature_names())]))
    structured = winner["factor_weight"] * factorized + (1.0 - winner["factor_weight"]) * direct
    observed_probabilities = winner["candidate_weight"] * structured + (1.0 - winner["candidate_weight"]) * baseline_observed
    predicted = np.asarray(LABELS)[np.argmax(observed_probabilities, axis=1)]
    rows = [
        {
            "indID1": pair[0], "indID2": pair[1], "family": family_indices[index],
            "relation": str(predicted[index]),
            "confidence": float(observed_probabilities[index].max()),
            "probabilities": {label: float(value) for label, value in zip(LABELS, observed_probabilities[index], strict=True)},
        }
        for index, pair in enumerate(observed_pairs)
    ]
    args.out_dir.mkdir(parents=True, exist_ok=True)
    joblib.dump({"bundle": bundle, "winner": winner}, args.out_dir / "model.joblib")
    (args.out_dir / "predictions.json").write_text(json.dumps(rows, indent=2) + "\n")
    summary = {
        "state": "DEVELOPMENT", "design": "degree-factorized histogram-gradient model with IBD-state and chromosome spectra",
        "selected": winner, "selection": screens,
        "development_evaluation": metrics(development, development_y),
        "development_baseline": metrics(baseline_development, development_y),
        "observed_pairs": len(rows), "observed_label_counts": dict(Counter(predicted)),
        "hbd_threshold_cm": args.hbd_threshold_cm, "submission_count": 0,
    }
    (args.out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({key: value for key, value in summary.items() if key not in {"selection"}}, indent=2))


if __name__ == "__main__":
    main()
