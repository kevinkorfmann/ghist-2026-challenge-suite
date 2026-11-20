from __future__ import annotations

import argparse
import copy
import csv
import json
from collections import Counter
from pathlib import Path

import joblib
import numpy as np
from scipy.stats import poisson

try:
    from scripts.fit_human_background_robust import (
        DEVELOPMENT_REGIMES,
        SELECTION_REGIMES,
        regime_datasets,
    )
    from scripts.fit_human_family_relationships import (
        BIN_EDGES,
        MEIOSES,
        ORDERED,
        background_correct,
        matrix,
        read_observed,
    )
    from scripts.validate_human_background_correction import (
        trim_shortest_background_shape,
    )
except ModuleNotFoundError:
    from fit_human_background_robust import (
        DEVELOPMENT_REGIMES,
        SELECTION_REGIMES,
        regime_datasets,
    )
    from fit_human_family_relationships import (
        BIN_EDGES,
        MEIOSES,
        ORDERED,
        background_correct,
        matrix,
        read_observed,
    )
    from validate_human_background_correction import trim_shortest_background_shape


def remove_shortest_segments(
    records: list[dict[str, object]],
    remove_count: int,
    minimum_cm: int,
    genome_cm: float,
) -> list[dict[str, object]]:
    output = copy.deepcopy(records)
    for row in output:
        lengths = [
            float(row[f"ordered_length_{rank + 1}"])
            for rank in range(ORDERED)
            if float(row[f"ordered_length_{rank + 1}"]) > 0
        ]
        eligible = sorted(value for value in lengths if value >= minimum_cm)
        removed = eligible[:remove_count]
        retained = list(lengths)
        for value in removed:
            retained.remove(value)
        retained.sort(reverse=True)
        for threshold in (5, 10, 20, 40):
            selected = [value for value in removed if value >= threshold]
            row[f"segments_ge_{threshold}"] = max(
                0.0, float(row[f"segments_ge_{threshold}"]) - len(selected)
            )
            total_key = f"total_cm_ge_{threshold}"
            row[total_key] = max(0.0, float(row[total_key]) - sum(selected))
            row[f"relatedness_ge_{threshold}"] = (
                float(row[total_key]) / (2.0 * genome_cm)
            )
        removed_bins, _ = np.histogram(np.asarray(removed), bins=np.asarray(BIN_EDGES))
        for index, count in enumerate(removed_bins):
            key = f"length_bin_{index}"
            row[key] = max(0.0, float(row[key]) - float(count))
        row["max_cm"] = retained[0] if retained else 0.0
        for rank in range(ORDERED):
            row[f"ordered_length_{rank + 1}"] = (
                retained[rank] if rank < len(retained) else 0.0
            )
    return output


def probability_stack(
    model,
    records: list[dict[str, object]],
    names: list[str],
    maximum_removed: int,
    minimum_cm: int,
    genome_cm: float,
) -> tuple[np.ndarray, np.ndarray]:
    probabilities = []
    for count in range(maximum_removed + 1):
        corrected = remove_shortest_segments(
            records, count, minimum_cm, genome_cm
        )
        probabilities.append(model.predict_proba(matrix(corrected, names)))
    return np.asarray(model.classes_), np.stack(probabilities, axis=0)


def combine_probabilities(
    stack: np.ndarray,
    mean_segments: float,
    confidence_power: float,
    mode: str,
) -> np.ndarray:
    counts = np.arange(stack.shape[0])
    prior = poisson.pmf(counts, mean_segments)
    prior[-1] += max(0.0, 1.0 - float(prior.sum()))
    confidence = np.max(stack, axis=2)
    weights = prior[:, None] * np.maximum(confidence, 1e-12) ** confidence_power
    if mode == "map":
        selected = np.argmax(weights, axis=0)
        return stack[selected, np.arange(stack.shape[1])]
    weights /= np.maximum(weights.sum(axis=0, keepdims=True), 1e-12)
    return np.sum(stack * weights[:, :, None], axis=0)


def prediction_metrics(
    classes: np.ndarray,
    probabilities: np.ndarray,
    records: list[dict[str, object]],
) -> dict[str, float]:
    truth = np.asarray([str(row["label"]) for row in records])
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


def current_baseline_metrics(
    model,
    names: list[str],
    records: list[dict[str, object]],
    reference: list[dict[str, object]],
    minimum_cm: int,
    genome_cm: float,
) -> dict[str, float]:
    corrected = copy.deepcopy(records)
    offsets = background_correct(corrected, reference, 0.25, genome_cm)
    trim_shortest_background_shape(
        corrected,
        float(offsets[f"segments_ge_{minimum_cm}"]),
        minimum_cm,
    )
    probabilities = model.predict_proba(matrix(corrected, names))
    return prediction_metrics(np.asarray(model.classes_), probabilities, records)


def candidate_specs(maximum_removed: int) -> list[dict[str, object]]:
    specs = [
        {"name": f"fixed_{count}", "mode": "fixed", "count": count}
        for count in range(maximum_removed + 1)
    ]
    for mean_segments in (1.5, 2.5, 3.5, 4.5, 5.5):
        for confidence_power in (0.0, 2.0, 5.0, 10.0):
            for mode in ("marginal", "map"):
                specs.append(
                    {
                        "name": f"{mode}_lambda_{mean_segments}_power_{confidence_power}",
                        "mode": mode,
                        "mean_segments": mean_segments,
                        "confidence_power": confidence_power,
                    }
                )
    return specs


def evaluate_specs(
    model,
    names: list[str],
    datasets,
    specs: list[dict[str, object]],
    maximum_removed: int,
    minimum_cm: int,
    genome_cm: float,
) -> list[dict[str, object]]:
    runs = {str(spec["name"]): [] for spec in specs}
    for regime, records in datasets:
        classes, stack = probability_stack(
            model,
            records,
            names,
            maximum_removed,
            minimum_cm,
            genome_cm,
        )
        for spec in specs:
            if spec["mode"] == "fixed":
                probabilities = stack[int(spec["count"])]
            else:
                probabilities = combine_probabilities(
                    stack,
                    float(spec["mean_segments"]),
                    float(spec["confidence_power"]),
                    str(spec["mode"]),
                )
            result = prediction_metrics(classes, probabilities, records)
            result["background_mean_segments"] = regime.mean_segments
            runs[str(spec["name"])].append(result)
    results = []
    for spec in specs:
        selected_runs = runs[str(spec["name"])]
        results.append(
            {
                **spec,
                "mean_exact_accuracy": float(
                    np.mean([row["exact_accuracy"] for row in selected_runs])
                ),
                "worst_exact_accuracy": float(
                    np.min([row["exact_accuracy"] for row in selected_runs])
                ),
                "mean_meioses_accuracy": float(
                    np.mean([row["meioses_accuracy"] for row in selected_runs])
                ),
                "runs": selected_runs,
            }
        )
    return results


def predict_observed(
    model,
    names: list[str],
    records: list[dict[str, object]],
    spec: dict[str, object],
    maximum_removed: int,
    minimum_cm: int,
    genome_cm: float,
) -> list[dict[str, object]]:
    classes, stack = probability_stack(
        model, records, names, maximum_removed, minimum_cm, genome_cm
    )
    if spec["mode"] == "fixed":
        probabilities = stack[int(spec["count"])]
    else:
        probabilities = combine_probabilities(
            stack,
            float(spec["mean_segments"]),
            float(spec["confidence_power"]),
            str(spec["mode"]),
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
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--selection", type=Path, required=True)
    parser.add_argument("--development-evaluation", type=Path, required=True)
    parser.add_argument("--observed", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--maximum-removed", type=int, default=8)
    parser.add_argument("--genome-cm", type=float, default=3346.298892538)
    parser.add_argument("--seed", type=int, default=122026)
    args = parser.parse_args()

    payload = joblib.load(args.model)
    model, names = payload["model"], payload["feature_names"]
    minimum_cm = min(
        int(name.rsplit("_", 1)[1])
        for name in names
        if name.startswith("segments_ge_")
    )
    selection = regime_datasets(
        args.selection, SELECTION_REGIMES, args.seed, args.genome_cm
    )
    specs = candidate_specs(args.maximum_removed)
    selection_results = evaluate_specs(
        model,
        names,
        selection,
        specs,
        args.maximum_removed,
        minimum_cm,
        args.genome_cm,
    )
    selection_results.sort(
        key=lambda row: (
            row["mean_exact_accuracy"],
            row["worst_exact_accuracy"],
            row["mean_meioses_accuracy"],
        ),
        reverse=True,
    )
    selected_spec = {
        key: value
        for key, value in selection_results[0].items()
        if key not in {"runs", "mean_exact_accuracy", "worst_exact_accuracy", "mean_meioses_accuracy"}
    }
    development = regime_datasets(
        args.development_evaluation,
        DEVELOPMENT_REGIMES,
        args.seed + 1000,
        args.genome_cm,
    )
    development_result = evaluate_specs(
        model,
        names,
        development,
        [selected_spec],
        args.maximum_removed,
        minimum_cm,
        args.genome_cm,
    )[0]
    observed = read_observed(args.observed)
    predictions = predict_observed(
        model,
        names,
        observed,
        selected_spec,
        args.maximum_removed,
        minimum_cm,
        args.genome_cm,
    )

    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "predictions.json").write_text(json.dumps(predictions, indent=2) + "\n")
    with (args.out_dir / "candidate.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("indID1", "indID2", "relation"))
        writer.writeheader()
        writer.writerows(
            {key: row[key] for key in writer.fieldnames} for row in predictions
        )
    summary = {
        "state": "DEVELOPMENT",
        "design": "ERSA-inspired marginal deconvolution of pair-specific background segments",
        "selected_spec": selected_spec,
        "selection": selection_results,
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
                "state": summary["state"],
                "selected_spec": selected_spec,
                "selection_winner": selection_results[0],
                "development_evaluation": development_result,
                "observed_label_counts": summary["observed_label_counts"],
                "submission_count": 0,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
