from __future__ import annotations

import argparse
import copy
import csv
import json
from collections import Counter
from pathlib import Path

import joblib
import numpy as np

try:
    from scripts.fit_human_background_robust import (
        DEVELOPMENT_REGIMES,
        SELECTION_REGIMES,
        augment_records,
        fit_direct_models,
        fit_hierarchical_model,
        predict_probabilities,
        random_regime,
        regime_datasets,
    )
    from scripts.fit_human_family_relationships import (
        LABELS,
        MEIOSES,
        background_correct,
        feature_names,
        matrix,
        read_observed,
    )
    from scripts.validate_human_background_correction import (
        read_pedsim_with_lengths,
        trim_shortest_background_shape,
    )
except ModuleNotFoundError:
    from fit_human_background_robust import (
        DEVELOPMENT_REGIMES,
        SELECTION_REGIMES,
        augment_records,
        fit_direct_models,
        fit_hierarchical_model,
        predict_probabilities,
        random_regime,
        regime_datasets,
    )
    from fit_human_family_relationships import (
        LABELS,
        MEIOSES,
        background_correct,
        feature_names,
        matrix,
        read_observed,
    )
    from validate_human_background_correction import (
        read_pedsim_with_lengths,
        trim_shortest_background_shape,
    )


def correct_records(
    records: list[dict[str, object]],
    reference: list[dict[str, object]],
    minimum_cm: int,
    genome_cm: float,
) -> list[dict[str, object]]:
    corrected = copy.deepcopy(records)
    offsets = background_correct(corrected, reference, 0.25, genome_cm)
    trim_shortest_background_shape(
        corrected,
        float(offsets[f"segments_ge_{minimum_cm}"]),
        minimum_cm,
    )
    return corrected


def corrected_training_set(
    paths: list[Path],
    augmentations: int,
    seed: int,
    minimum_cm: int,
    genome_cm: float,
) -> tuple[list[dict[str, object]], list[dict[str, float]]]:
    rng = np.random.default_rng(seed)
    output = []
    regimes = []
    for path in paths:
        reference, lengths = read_pedsim_with_lengths(path, genome_cm)
        output.extend(reference)
        for _ in range(augmentations):
            regime = random_regime(rng)
            augmented = augment_records(reference, lengths, regime, rng, genome_cm)
            output.extend(correct_records(augmented, reference, minimum_cm, genome_cm))
            regimes.append(
                {
                    "mean_segments": regime.mean_segments,
                    "excess_scale_cm": regime.excess_scale_cm,
                    "pair_heterogeneity": regime.pair_heterogeneity,
                    "long_tail_probability": regime.long_tail_probability,
                    "long_tail_scale": regime.long_tail_scale,
                }
            )
    return output, regimes


def corrected_regime_datasets(
    path: Path,
    regimes,
    seed: int,
    minimum_cm: int,
    genome_cm: float,
):
    reference, _ = read_pedsim_with_lengths(path, genome_cm)
    return [
        (
            regime,
            correct_records(records, reference, minimum_cm, genome_cm),
        )
        for regime, records in regime_datasets(path, regimes, seed, genome_cm)
    ]


def align_probabilities(
    classes: np.ndarray, probabilities: np.ndarray, target: np.ndarray
) -> np.ndarray:
    output = np.zeros((probabilities.shape[0], len(target)), dtype=float)
    index = {str(label): offset for offset, label in enumerate(target)}
    for offset, label in enumerate(classes):
        output[:, index[str(label)]] = probabilities[:, offset]
    return output


def metrics(
    classes: np.ndarray,
    probabilities: np.ndarray,
    records: list[dict[str, object]],
) -> dict[str, float]:
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


def evaluate_blend(
    corrected_payload,
    baseline_model,
    baseline_names,
    datasets,
    corrected_names,
    weight: float,
) -> dict[str, object]:
    classes = np.asarray(LABELS)
    runs = []
    for regime, records in datasets:
        corrected_classes, corrected_probabilities = predict_probabilities(
            corrected_payload, matrix(records, corrected_names)
        )
        baseline_probabilities = baseline_model.predict_proba(
            matrix(records, baseline_names)
        )
        corrected_aligned = align_probabilities(
            corrected_classes, corrected_probabilities, classes
        )
        baseline_aligned = align_probabilities(
            np.asarray(baseline_model.classes_), baseline_probabilities, classes
        )
        probabilities = (
            weight * corrected_aligned + (1.0 - weight) * baseline_aligned
        )
        result = metrics(classes, probabilities, records)
        result["background_mean_segments"] = regime.mean_segments
        runs.append(result)
    return {
        "mean_exact_accuracy": float(np.mean([row["exact_accuracy"] for row in runs])),
        "worst_exact_accuracy": float(np.min([row["exact_accuracy"] for row in runs])),
        "mean_meioses_accuracy": float(np.mean([row["meioses_accuracy"] for row in runs])),
        "runs": runs,
    }


def observed_predictions(
    corrected_payload,
    baseline_model,
    baseline_names,
    records,
    corrected_names,
    weight: float,
):
    classes = np.asarray(LABELS)
    corrected_classes, corrected_probabilities = predict_probabilities(
        corrected_payload, matrix(records, corrected_names)
    )
    baseline_probabilities = baseline_model.predict_proba(matrix(records, baseline_names))
    probabilities = weight * align_probabilities(
        corrected_classes, corrected_probabilities, classes
    ) + (1.0 - weight) * align_probabilities(
        np.asarray(baseline_model.classes_), baseline_probabilities, classes
    )
    predicted = classes[np.argmax(probabilities, axis=1)]
    return [
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
        for index, record in enumerate(records)
    ]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline-model", type=Path, required=True)
    parser.add_argument("--calibration", type=Path, nargs="+", required=True)
    parser.add_argument("--selection", type=Path, required=True)
    parser.add_argument("--development-evaluation", type=Path, required=True)
    parser.add_argument("--observed", type=Path, required=True)
    parser.add_argument("--background-reference", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--minimum-cm", type=int, default=20)
    parser.add_argument("--augmentations", type=int, default=4)
    parser.add_argument("--genome-cm", type=float, default=3346.298892538)
    parser.add_argument("--seed", type=int, default=122126)
    args = parser.parse_args()

    baseline = joblib.load(args.baseline_model)
    names = feature_names(args.minimum_cm)
    training, training_regimes = corrected_training_set(
        args.calibration,
        args.augmentations,
        args.seed,
        args.minimum_cm,
        args.genome_cm,
    )
    train_x = matrix(training, names)
    train_y = np.asarray([str(record["label"]) for record in training])
    candidates = fit_direct_models(train_x, train_y, args.seed)
    candidates.append(fit_hierarchical_model(train_x, train_y, args.seed + 100))
    selection = corrected_regime_datasets(
        args.selection,
        SELECTION_REGIMES,
        args.seed + 1000,
        args.minimum_cm,
        args.genome_cm,
    )
    results = []
    fitted = {}
    for name, payload in candidates:
        fitted[name] = payload
        for weight in (0.25, 0.50, 0.75, 1.0):
            result = evaluate_blend(
                payload,
                baseline["model"],
                baseline["feature_names"],
                selection,
                names,
                weight,
            )
            result.update({"name": name, "corrected_model_weight": weight})
            results.append(result)
    results.sort(
        key=lambda row: (
            row["mean_exact_accuracy"],
            row["worst_exact_accuracy"],
            row["mean_meioses_accuracy"],
        ),
        reverse=True,
    )
    winner = results[0]
    selected_payload = fitted[str(winner["name"])]
    weight = float(winner["corrected_model_weight"])
    development = corrected_regime_datasets(
        args.development_evaluation,
        DEVELOPMENT_REGIMES,
        args.seed + 2000,
        args.minimum_cm,
        args.genome_cm,
    )
    development_result = evaluate_blend(
        selected_payload,
        baseline["model"],
        baseline["feature_names"],
        development,
        names,
        weight,
    )
    observed = read_observed(args.observed)
    reference = read_observed(args.background_reference)
    corrected_observed = correct_records(
        observed, reference, args.minimum_cm, args.genome_cm
    )
    predictions = observed_predictions(
        selected_payload,
        baseline["model"],
        baseline["feature_names"],
        corrected_observed,
        names,
        weight,
    )

    args.out_dir.mkdir(parents=True, exist_ok=True)
    joblib.dump(
        {
            "corrected_payload": selected_payload,
            "feature_names": names,
            "corrected_model_weight": weight,
        },
        args.out_dir / "model.joblib",
    )
    (args.out_dir / "predictions.json").write_text(json.dumps(predictions, indent=2) + "\n")
    with (args.out_dir / "candidate.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("indID1", "indID2", "relation"))
        writer.writeheader()
        writer.writerows(
            {key: row[key] for key in writer.fieldnames} for row in predictions
        )
    summary = {
        "state": "DEVELOPMENT",
        "design": "domain-adapted learning on globally corrected background contamination with clean-model blending",
        "training_pairs": len(training),
        "training_regimes": training_regimes,
        "selection": results,
        "selected_model": winner["name"],
        "corrected_model_weight": weight,
        "development_evaluation": development_result,
        "observed_pairs": len(predictions),
        "observed_label_counts": dict(Counter(row["relation"] for row in predictions)),
        "final_holdout_evaluated": False,
        "submission_count": 0,
    }
    (args.out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(
        json.dumps(
            {
                "selected_model": summary["selected_model"],
                "corrected_model_weight": weight,
                "selection_winner": winner,
                "development_evaluation": development_result,
                "observed_label_counts": summary["observed_label_counts"],
                "submission_count": 0,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
