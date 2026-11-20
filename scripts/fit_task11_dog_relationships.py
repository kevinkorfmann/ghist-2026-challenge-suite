from __future__ import annotations

import argparse
import gzip
import json
from collections import defaultdict
from pathlib import Path

from ghist_relatedness.validation import validate_simulation_splits

import joblib
import numpy as np
from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier
from sklearn.metrics import accuracy_score, balanced_accuracy_score, confusion_matrix

try:
    from scripts.reconstruct_bonsai_pedigree import convert_segments, read_map
except ModuleNotFoundError:
    from reconstruct_bonsai_pedigree import convert_segments, read_map


LABELS = (
    "siblings",
    "half_sibling",
    "avuncular",
    "half_avuncular",
    "avuncular_once_removed",
    "avuncular_twice_removed",
    "full_1st_cousin",
    "half_1st_cousin",
)
THRESHOLDS = (2.0, 5.0, 10.0, 20.0, 40.0)
TOP_LENGTHS = 12


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


def feature_names() -> list[str]:
    names = [
        "segments",
        "chromosomes",
        "total_ibd1",
        "total_ibd2",
        "total_sharing",
        "ibd2_fraction",
        "mean_length",
        "median_length",
        "q75_length",
        "q90_length",
        "max_length",
    ]
    for threshold in THRESHOLDS:
        tag = str(threshold).replace(".", "p")
        names.extend(
            (
                f"segments_ge_{tag}",
                f"ibd1_ge_{tag}",
                f"ibd2_ge_{tag}",
                f"sharing_ge_{tag}",
            )
        )
    names.extend(f"ordered_length_{rank}" for rank in range(1, TOP_LENGTHS + 1))
    return names


def summarize_segments(
    records: list[tuple[str, float, bool]], exposure_scale: float = 1.0
) -> np.ndarray:
    lengths = np.asarray([record[1] for record in records], dtype=float)
    full = np.asarray([record[2] for record in records], dtype=bool)
    chromosomes = {record[0] for record in records}
    if lengths.size == 0:
        return np.zeros(len(feature_names()), dtype=float)
    total_ibd1 = float(lengths[~full].sum()) * exposure_scale
    total_ibd2 = float(lengths[full].sum()) * exposure_scale
    total_sharing = total_ibd1 + 2.0 * total_ibd2
    values = [
        float(lengths.size),
        float(len(chromosomes)),
        total_ibd1,
        total_ibd2,
        total_sharing,
        2.0 * total_ibd2 / max(total_sharing, 1e-12),
        float(lengths.mean()),
        float(np.median(lengths)),
        float(np.quantile(lengths, 0.75)),
        float(np.quantile(lengths, 0.90)),
        float(lengths.max()),
    ]
    for threshold in THRESHOLDS:
        selected = lengths >= threshold
        ibd1 = float(lengths[selected & ~full].sum()) * exposure_scale
        ibd2 = float(lengths[selected & full].sum()) * exposure_scale
        values.extend((float(selected.sum()), ibd1, ibd2, ibd1 + 2.0 * ibd2))
    ordered = np.sort(lengths)[::-1][:TOP_LENGTHS]
    values.extend(np.pad(ordered, (0, TOP_LENGTHS - ordered.size)).tolist())
    return np.asarray(values, dtype=float)


def read_pedsim(path: Path) -> tuple[np.ndarray, np.ndarray]:
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
        np.vstack([summarize_segments(grouped[pair]) for pair in pairs]),
        np.asarray([labels[pair] for pair in pairs]),
    )


def read_clusters(path: Path) -> list[list[str]]:
    payload = json.loads(path.read_text())
    families = payload.get("families", payload.get("clusters"))
    if isinstance(families, dict):
        families = list(families.values())
    if not isinstance(families, list):
        raise ValueError("Family cluster file has no families")
    return [sorted(family) for family in families]


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    merged: list[tuple[int, int]] = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def subtract_intervals(
    start: int, end: int, masks: list[tuple[int, int]]
) -> list[tuple[int, int]]:
    pieces = [(start, end)]
    for mask_start, mask_end in masks:
        updated = []
        for left, right in pieces:
            if mask_end <= left or mask_start >= right:
                updated.append((left, right))
            else:
                if left < mask_start:
                    updated.append((left, mask_start))
                if mask_end < right:
                    updated.append((mask_end, right))
        pieces = updated
    return pieces


def read_hbd_masks(
    path: Path, threshold_cm: float
) -> dict[tuple[str, str], list[tuple[int, int]]]:
    raw: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
    with gzip.open(path, "rt") as handle:
        for line in handle:
            sample, _, _, _, chromosome, start, end, length = line.split()
            if float(length) >= threshold_cm:
                raw[(sample, chromosome)].append((int(start), int(end)))
    return {key: merge_intervals(intervals) for key, intervals in raw.items()}


def read_observed(
    ibd_path: Path,
    hbd_path: Path,
    map_path: Path,
    cluster_path: Path,
    hbd_threshold_cm: float,
    summary_function=summarize_segments,
) -> tuple[np.ndarray, list[tuple[str, str]], list[int]]:
    families = read_clusters(cluster_path)
    sample_to_family = {
        sample: family_index
        for family_index, family in enumerate(families)
        for sample in family
    }
    allowed = set(sample_to_family)
    segments, _, id_to_sample = convert_segments(
        ibd_path, map_path, threshold_cm=0.0, allowed_samples=allowed
    )
    genetic_map = read_map(map_path)
    genome_cm = float(sum(values[1][-1] - values[1][0] for values in genetic_map.values()))
    hbd_masks = read_hbd_masks(hbd_path, hbd_threshold_cm)
    grouped: dict[tuple[str, str], list[tuple[str, float, bool]]] = defaultdict(list)
    pair_masks: dict[tuple[str, str, str], list[tuple[int, int]]] = {}
    pair_scales: dict[tuple[str, str], float] = {}
    for left, right, chromosome, start, end, is_full, _ in segments:
        pair = tuple(sorted((id_to_sample[int(left)], id_to_sample[int(right)])))
        if sample_to_family[pair[0]] == sample_to_family[pair[1]]:
            mask_key = (pair[0], pair[1], str(chromosome))
            if mask_key not in pair_masks:
                pair_masks[mask_key] = merge_intervals(
                    hbd_masks.get((pair[0], str(chromosome)), [])
                    + hbd_masks.get((pair[1], str(chromosome)), [])
                )
            positions, cm = genetic_map[str(chromosome)]
            for piece_start, piece_end in subtract_intervals(
                int(start), int(end), pair_masks[mask_key]
            ):
                bounds = np.interp(
                    np.asarray([piece_start, piece_end], dtype=float), positions, cm
                )
                length = float(bounds[1] - bounds[0])
                if length > 0:
                    grouped[pair].append((str(chromosome), length, bool(is_full)))
    pairs = []
    family_indices = []
    for family_index, family in enumerate(families):
        for offset, left in enumerate(family):
            for right in family[offset + 1 :]:
                pairs.append((left, right))
                family_indices.append(family_index)
    for pair in pairs:
        masked_cm = 0.0
        for chromosome, (positions, cm) in genetic_map.items():
            masks = merge_intervals(
                hbd_masks.get((pair[0], chromosome), [])
                + hbd_masks.get((pair[1], chromosome), [])
            )
            for start, end in masks:
                bounds = np.interp(np.asarray([start, end], dtype=float), positions, cm)
                masked_cm += float(bounds[1] - bounds[0])
        pair_scales[pair] = genome_cm / max(genome_cm - masked_cm, 0.25 * genome_cm)
    matrix = np.vstack(
        [
            summary_function(grouped[pair], exposure_scale=pair_scales[pair])
            for pair in pairs
        ]
    )
    return matrix, pairs, family_indices


def candidate_models(seed: int):
    for leaf in (1, 2, 4):
        yield (
            f"extra_trees_leaf_{leaf}",
            ExtraTreesClassifier(
                n_estimators=800,
                min_samples_leaf=leaf,
                max_features=0.75,
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
                max_features=0.75,
                class_weight="balanced_subsample",
                n_jobs=-1,
                random_state=seed,
            ),
        )


def score(model, matrix: np.ndarray, truth: np.ndarray) -> dict[str, object]:
    prediction = model.predict(matrix)
    return {
        "accuracy": float(accuracy_score(truth, prediction)),
        "balanced_accuracy": float(balanced_accuracy_score(truth, prediction)),
        "confusion": confusion_matrix(truth, prediction, labels=LABELS).tolist(),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--calibration", type=Path, nargs="+", required=True)
    parser.add_argument("--validation", type=Path, required=True)
    parser.add_argument("--evaluation", type=Path, required=True)
    parser.add_argument("--observed-ibd", type=Path, required=True)
    parser.add_argument("--observed-hbd", type=Path, required=True)
    parser.add_argument("--hbd-threshold-cm", type=float, default=5.0)
    parser.add_argument("--map", type=Path, required=True)
    parser.add_argument("--clusters", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=260824)
    args = parser.parse_args()
    validate_simulation_splits({
        'calibration': args.calibration,
        'validation': [args.validation],
        'evaluation': [args.evaluation],
    })

    training = [read_pedsim(path) for path in args.calibration]
    train_x = np.vstack([item[0] for item in training])
    train_y = np.concatenate([item[1] for item in training])
    validation_x, validation_y = read_pedsim(args.validation)
    evaluation_x, evaluation_y = read_pedsim(args.evaluation)

    candidates = []
    fitted = []
    for name, model in candidate_models(args.seed):
        model.fit(train_x, train_y)
        result = score(model, validation_x, validation_y)
        result["name"] = name
        candidates.append(result)
        fitted.append((name, model))
    candidates.sort(key=lambda item: (item["balanced_accuracy"], item["accuracy"]), reverse=True)
    selected_name = candidates[0]["name"]
    selected = dict(fitted)[selected_name]
    evaluation = score(selected, evaluation_x, evaluation_y)

    observed_x, observed_pairs, family_indices = read_observed(
        args.observed_ibd,
        args.observed_hbd,
        args.map,
        args.clusters,
        args.hbd_threshold_cm,
    )
    probabilities = selected.predict_proba(observed_x)
    predicted = selected.classes_[np.argmax(probabilities, axis=1)]
    rows = []
    for index, pair in enumerate(observed_pairs):
        order = np.argsort(probabilities[index])[::-1]
        rows.append(
            {
                "indID1": pair[0],
                "indID2": pair[1],
                "family": family_indices[index],
                "relation": str(predicted[index]),
                "confidence": float(probabilities[index, order[0]]),
                "runner_up": str(selected.classes_[order[1]]),
                "runner_up_probability": float(probabilities[index, order[1]]),
                "probabilities": {
                    str(label): float(probability)
                    for label, probability in zip(
                        selected.classes_, probabilities[index], strict=True
                    )
                },
                "features": {
                    name: float(value)
                    for name, value in zip(feature_names(), observed_x[index])
                },
            }
        )

    args.out_dir.mkdir(parents=True, exist_ok=True)
    joblib.dump(selected, args.out_dir / "model.joblib")
    (args.out_dir / "predictions.json").write_text(json.dumps(rows, indent=2) + "\n")
    summary = {
        "selected_model": selected_name,
        "feature_names": feature_names(),
        "labels": list(LABELS),
        "training_pairs": int(train_y.size),
        "validation_pairs": int(validation_y.size),
        "evaluation_pairs": int(evaluation_y.size),
        "model_selection": candidates,
        "fresh_evaluation": evaluation,
        "observed_pairs": len(rows),
        "hbd_threshold_cm": args.hbd_threshold_cm,
        "observed_label_counts": {
            label: sum(row["relation"] == label for row in rows) for label in LABELS
        },
        "submission_count": 0,
    }
    (args.out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
