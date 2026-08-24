from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import joblib
import numpy as np
from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier
from sklearn.metrics import accuracy_score, balanced_accuracy_score, confusion_matrix

from ghist_relatedness.labels import HUMAN_LABELS
try:
    from scripts.summarize_hapibd_shape import BIN_EDGES, ORDERED
except ModuleNotFoundError:
    from summarize_hapibd_shape import BIN_EDGES, ORDERED


LABELS = tuple(dict.fromkeys(HUMAN_LABELS.values()))
MEIOSES = {label: up + down for (up, down, _), label in HUMAN_LABELS.items()}
THRESHOLDS = (5, 10, 20, 40)


def relationship_label(sample: str) -> str:
    prefix = sample.split("_g", 1)[0]
    matches = [
        label
        for label in LABELS
        if prefix.startswith(label) and prefix[len(label) :].isdigit()
    ]
    if len(matches) != 1:
        raise ValueError(f"Cannot recover relationship label from {sample}")
    return matches[0]


def summarize_raw_lengths(
    records: list[tuple[str, float]], genome_cm: float
) -> dict[str, float]:
    lengths = np.asarray([record[1] for record in records], dtype=float)
    output = {
        "max_cm": float(lengths.max()) if lengths.size else 0.0,
    }
    for threshold in THRESHOLDS:
        selected = lengths[lengths >= threshold]
        total = float(selected.sum())
        output[f"segments_ge_{threshold}"] = float(selected.size)
        output[f"total_cm_ge_{threshold}"] = total
        output[f"relatedness_ge_{threshold}"] = total / (2.0 * genome_cm)
    counts, _ = np.histogram(lengths, bins=np.asarray(BIN_EDGES))
    output.update(
        {f"length_bin_{index}": float(count) for index, count in enumerate(counts)}
    )
    ordered = np.sort(lengths)[::-1][:ORDERED]
    output.update(
        {
            f"ordered_length_{rank + 1}": float(ordered[rank])
            if rank < ordered.size
            else 0.0
            for rank in range(ORDERED)
        }
    )
    return output


def read_pedsim(path: Path, genome_cm: float) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str], list[tuple[str, float]]] = defaultdict(list)
    labels: dict[tuple[str, str], str] = {}
    with path.open() as handle:
        for line in handle:
            left, right, chromosome, _, _, state, _, _, length = line.split()
            pair = tuple(sorted((left, right)))
            label = relationship_label(left)
            if relationship_label(right) != label:
                raise ValueError(f"Mixed relationship labels for {pair}")
            labels[pair] = label
            grouped[pair].append((chromosome, float(length)))
            if state == "IBD2":
                grouped[pair].append((chromosome, float(length)))
    output = []
    for pair in sorted(grouped):
        record: dict[str, object] = {
            "indID1": pair[0],
            "indID2": pair[1],
            "label": labels[pair],
        }
        record.update(summarize_raw_lengths(grouped[pair], genome_cm))
        output.append(record)
    return output


def read_observed(path: Path) -> list[dict[str, object]]:
    with path.open() as handle:
        rows = list(csv.DictReader(handle))
    return [
        {
            key: value if key in {"indID1", "indID2"} else float(value)
            for key, value in row.items()
        }
        for row in rows
    ]


def background_correct(
    observed: list[dict[str, object]],
    reference: list[dict[str, object]],
    quantile: float,
    genome_cm: float,
) -> dict[str, float]:
    offsets = {}
    for threshold in THRESHOLDS:
        count_key = f"segments_ge_{threshold}"
        rel_key = f"relatedness_ge_{threshold}"
        count_offset = max(
            0.0,
            float(np.quantile([row[count_key] for row in observed], quantile))
            - float(np.quantile([row[count_key] for row in reference], quantile)),
        )
        rel_offset = max(
            0.0,
            float(np.quantile([row[rel_key] for row in observed], quantile))
            - float(np.quantile([row[rel_key] for row in reference], quantile)),
        )
        offsets[count_key] = count_offset
        offsets[rel_key] = rel_offset
        offsets[f"total_cm_ge_{threshold}"] = rel_offset * 2.0 * genome_cm
    for row in observed:
        for key, offset in offsets.items():
            row[key] = max(0.0, float(row[key]) - offset)
    return offsets


def feature_names(minimum_cm: int) -> list[str]:
    names = ["max_cm"]
    for threshold in THRESHOLDS:
        if threshold >= minimum_cm:
            names.extend(
                (
                    f"segments_ge_{threshold}",
                    f"total_cm_ge_{threshold}",
                    f"relatedness_ge_{threshold}",
                )
            )
    for index, lower in enumerate(BIN_EDGES[:-1]):
        if lower >= minimum_cm:
            names.append(f"length_bin_{index}")
    names.extend(f"ordered_length_{rank + 1}" for rank in range(ORDERED))
    return names


def matrix(records: list[dict[str, object]], names: list[str]) -> np.ndarray:
    minimums = [
        int(name.rsplit("_", 1)[1])
        for name in names
        if name.startswith("segments_ge_")
    ]
    minimum = min(minimums) if minimums else 40
    values = []
    for record in records:
        row = []
        for name in names:
            value = float(record[name])
            if name.startswith("ordered_length_") and value < minimum:
                value = 0.0
            row.append(value)
        values.append(row)
    return np.asarray(values)


def candidate_models(seed: int):
    for leaf in (2, 4):
        yield (
            f"extra_trees_leaf_{leaf}",
            ExtraTreesClassifier(
                n_estimators=800,
                min_samples_leaf=leaf,
                max_features=0.8,
                class_weight="balanced",
                n_jobs=-1,
                random_state=seed,
            ),
        )
        yield (
            f"random_forest_leaf_{leaf}",
            RandomForestClassifier(
                n_estimators=800,
                min_samples_leaf=leaf,
                max_features=0.8,
                class_weight="balanced_subsample",
                n_jobs=-1,
                random_state=seed,
            ),
        )


def score(model, records: list[dict[str, object]], names: list[str]) -> dict[str, object]:
    truth = np.asarray([record["label"] for record in records])
    predicted = model.predict(matrix(records, names))
    truth_meioses = np.asarray([MEIOSES[str(label)] for label in truth])
    predicted_meioses = np.asarray([MEIOSES[str(label)] for label in predicted])
    return {
        "accuracy": float(accuracy_score(truth, predicted)),
        "balanced_accuracy": float(balanced_accuracy_score(truth, predicted)),
        "meioses_accuracy": float(accuracy_score(truth_meioses, predicted_meioses)),
        "confusion": confusion_matrix(truth, predicted, labels=LABELS).tolist(),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--calibration", type=Path, nargs="+", required=True)
    parser.add_argument("--validation", type=Path, required=True)
    parser.add_argument("--evaluation", type=Path, required=True)
    parser.add_argument("--observed", type=Path, required=True)
    parser.add_argument("--background-reference", type=Path)
    parser.add_argument("--background-quantile", type=float, default=0.5)
    parser.add_argument("--minimum-cm", type=int, choices=THRESHOLDS, default=5)
    parser.add_argument("--genome-cm", type=float, default=3346.298892538)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=101226)
    args = parser.parse_args()

    training = [read_pedsim(path, args.genome_cm) for path in args.calibration]
    train = [record for replicate in training for record in replicate]
    validation = read_pedsim(args.validation, args.genome_cm)
    evaluation = read_pedsim(args.evaluation, args.genome_cm)
    observed = read_observed(args.observed)
    offsets = {}
    if args.background_reference is not None:
        offsets = background_correct(
            observed,
            read_observed(args.background_reference),
            args.background_quantile,
            args.genome_cm,
        )
    names = feature_names(args.minimum_cm)
    train_x = matrix(train, names)
    train_y = np.asarray([record["label"] for record in train])

    candidates = []
    fitted = []
    for name, model in candidate_models(args.seed):
        model.fit(train_x, train_y)
        result = score(model, validation, names)
        result["name"] = name
        candidates.append(result)
        fitted.append((name, model))
    candidates.sort(
        key=lambda item: (
            item["meioses_accuracy"],
            item["balanced_accuracy"],
            item["accuracy"],
        ),
        reverse=True,
    )
    selected_name = candidates[0]["name"]
    selected = dict(fitted)[selected_name]
    evaluation_score = score(selected, evaluation, names)

    probabilities = selected.predict_proba(matrix(observed, names))
    predictions = selected.classes_[np.argmax(probabilities, axis=1)]
    rows = []
    for index, record in enumerate(observed):
        rows.append(
            {
                "indID1": record["indID1"],
                "indID2": record["indID2"],
                "relation": str(predictions[index]),
                "meioses_count": MEIOSES[str(predictions[index])],
                "confidence": float(np.max(probabilities[index])),
                "probabilities": {
                    str(label): float(probability)
                    for label, probability in zip(
                        selected.classes_, probabilities[index], strict=True
                    )
                },
            }
        )

    args.out_dir.mkdir(parents=True, exist_ok=True)
    joblib.dump(
        {"model": selected, "feature_names": names}, args.out_dir / "model.joblib"
    )
    (args.out_dir / "predictions.json").write_text(json.dumps(rows, indent=2) + "\n")
    summary = {
        "selected_model": selected_name,
        "minimum_cm": args.minimum_cm,
        "background_quantile": args.background_quantile
        if args.background_reference is not None
        else None,
        "background_offsets": offsets,
        "feature_names": names,
        "labels": list(LABELS),
        "training_pairs": len(train),
        "validation_pairs": len(validation),
        "evaluation_pairs": len(evaluation),
        "model_selection": candidates,
        "fresh_evaluation": evaluation_score,
        "observed_pairs": len(rows),
        "observed_label_counts": {
            label: sum(row["relation"] == label for row in rows) for label in LABELS
        },
        "submission_count": 0,
    }
    (args.out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
