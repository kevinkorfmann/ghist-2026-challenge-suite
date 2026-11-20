"""Independent descriptive scoring of complete held-out relatedness cohorts.

Missing or abstained pairs remain in every accuracy denominator. Bootstrap
units are whole cohorts, never individual related pairs. This tool does not
select models, establish cohort independence, or authorize a release.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from numbers import Real
from pathlib import Path

import numpy as np


def _pair(row: dict) -> tuple[str, str]:
    ids = (row["indID1"], row["indID2"])
    if any(not isinstance(value, str) or not value.strip() for value in ids):
        raise ValueError("Pair identifiers must be nonempty strings")
    if ids[0] == ids[1]:
        raise ValueError("Self pairs are not allowed")
    return tuple(sorted(ids))


def _index(rows: list[dict], labels: list[str], *, truth: bool) -> dict:
    output = {}
    for row in rows:
        pair = _pair(row)
        if pair in output:
            raise ValueError(f"Duplicate unordered pair: {pair}")
        relation = row["relation"]
        if relation not in labels and (truth or relation is not None):
            raise ValueError(f"Unsupported relationship: {relation!r}")
        probabilities = row.get("probabilities")
        if probabilities is not None:
            if not isinstance(probabilities, dict) or set(probabilities) != set(labels):
                raise ValueError("Probabilities must cover exactly the label vocabulary")
            if any(
                isinstance(value, bool)
                or not isinstance(value, Real)
                or not math.isfinite(value)
                or not 0 <= value <= 1
                for value in probabilities.values()
            ):
                raise ValueError("Invalid probability value")
            if not math.isclose(sum(probabilities.values()), 1.0, abs_tol=1e-8):
                raise ValueError("Probabilities must sum to one")
        output[pair] = row
    return output


def _metrics(confusion: np.ndarray, labels: list[str]) -> dict:
    denominators = confusion.sum(axis=1)
    correct = np.diag(confusion[:, :-1])
    recalls = np.divide(
        correct, denominators,
        out=np.full(len(labels), np.nan), where=denominators > 0,
    )
    represented = denominators > 0
    return {
        "pairs": int(denominators.sum()),
        "correct": int(correct.sum()),
        "exact_accuracy": float(correct.sum() / denominators.sum()),
        "balanced_accuracy": float(recalls.mean()) if represented.all() else None,
        "observed_label_macro_recall": float(recalls[represented].mean()),
        "all_labels_represented": bool(represented.all()),
        "per_label": {
            label: {
                "pairs": int(denominators[i]),
                "correct": int(correct[i]),
                "recall": float(recalls[i]) if represented[i] else None,
            }
            for i, label in enumerate(labels)
        },
        "confusion": confusion.tolist(),
        "confusion_rows": labels,
        "confusion_columns": labels + ["missing_or_abstained"],
    }


def _calibration(confidences: list[float], correct: list[bool]) -> dict:
    confidence = np.asarray(confidences)
    correctness = np.asarray(correct, dtype=float)
    bins = []
    ece = 0.0
    # Fixed 10 equal-width bins: [0,.1), ..., [.9,1].
    indices = np.minimum((confidence * 10).astype(int), 9)
    for index in range(10):
        selected = indices == index
        count = int(selected.sum())
        mean_confidence = float(confidence[selected].mean()) if count else None
        accuracy = float(correctness[selected].mean()) if count else None
        if count:
            ece += count / len(confidence) * abs(mean_confidence - accuracy)
        bins.append({
            "lower": index / 10, "upper": (index + 1) / 10,
            "pairs": count, "mean_confidence": mean_confidence, "accuracy": accuracy,
        })
    return {
        "ece": ece,
        "bins": bins,
        "pairs": len(confidence),
        "confidence_definition": "Reported probability of the actual emitted label",
        "scope": "Descriptive pair-weighted calibration; ECE alone is not a release gate",
    }


def evaluate_cohorts(
    cohorts: list[dict], labels: list[str], *,
    bootstrap_replicates: int = 2000, bootstrap_seed: int = 20261120,
) -> dict:
    if (
        not labels
        or any(not isinstance(label, str) or not label for label in labels)
        or len(set(labels)) != len(labels)
    ):
        raise ValueError("Supply a nonempty unique label vocabulary")
    if not cohorts:
        raise ValueError("At least one cohort is required")
    if (
        isinstance(bootstrap_replicates, bool)
        or not isinstance(bootstrap_replicates, int)
        or bootstrap_replicates < 1
    ):
        raise ValueError("bootstrap_replicates must be a positive integer")
    identifiers = set()
    matrices = []
    summaries = []
    confidence = []
    correctness = []
    complete_probability_coverage = True
    for cohort in cohorts:
        identifier = cohort["cohort_id"]
        if not isinstance(identifier, str) or not identifier or identifier in identifiers:
            raise ValueError("Cohort identifiers must be nonempty and unique")
        identifiers.add(identifier)
        truth = _index(cohort["truth"], labels, truth=True)
        predictions = _index(cohort["predictions"], labels, truth=False)
        if not truth:
            raise ValueError("A cohort must contain truth pairs")
        missing = len(truth.keys() - predictions.keys())
        extras = len(predictions.keys() - truth.keys())
        abstentions = 0
        matrix = np.zeros((len(labels), len(labels) + 1), dtype=np.int64)
        for pair, expected in truth.items():
            predicted = predictions.get(pair)
            relation = predicted["relation"] if predicted else None
            abstentions += int(predicted is not None and relation is None)
            column = labels.index(relation) if relation is not None else len(labels)
            matrix[labels.index(expected["relation"]), column] += 1
            if relation is None or predicted.get("probabilities") is None:
                complete_probability_coverage = False
            else:
                confidence.append(predicted["probabilities"][relation])
                correctness.append(relation == expected["relation"])
        complete_probability_coverage &= extras == 0
        summary = _metrics(matrix, labels)
        summary.update({
            "cohort_id": identifier, "missing_pairs": missing,
            "extra_pairs": extras, "abstained_pairs": abstentions,
            "pair_coverage_passed": missing == extras == abstentions == 0,
        })
        summaries.append(summary)
        matrices.append(matrix)
    stacked = np.stack(matrices)
    pooled = _metrics(stacked.sum(axis=0), labels)
    rng = np.random.default_rng(bootstrap_seed)
    bootstrap_values = []
    omitted = 0
    if len(cohorts) >= 2:
        for _ in range(bootstrap_replicates):
            indices = rng.integers(0, len(cohorts), size=len(cohorts))
            resampled = stacked[indices].sum(axis=0)
            denominators = resampled.sum(axis=1)
            if np.any(denominators == 0):
                omitted += 1
                continue
            correct = np.diag(resampled[:, :-1])
            bootstrap_values.append((
                correct.sum() / denominators.sum(), (correct / denominators).mean()
            ))
    bootstrap = {
        "unit": "whole cohort",
        "seed": bootstrap_seed,
        "replicates_requested": bootstrap_replicates,
        "replicates_evaluated": len(bootstrap_values),
        "replicates_missing_labels": omitted,
        "confidence_level": 0.95,
        "exact_accuracy_interval": None,
        "balanced_accuracy_interval": None,
        "status": "insufficient_cohorts" if len(cohorts) < 2 else "missing_label_support",
        "limitations": "Percentile bootstrap is descriptive; few cohorts or unrepresentative models cannot establish generalization.",
    }
    # Do not silently condition intervals on bootstrap samples retaining all labels.
    if len(cohorts) >= 2 and not omitted:
        bounds = np.quantile(np.asarray(bootstrap_values), [0.025, 0.975], axis=0)
        bootstrap.update({
            "status": "computed",
            "exact_accuracy_interval": bounds[:, 0].tolist(),
            "balanced_accuracy_interval": bounds[:, 1].tolist(),
        })
    return {
        "schema_version": "relatedness-cohort-score-v1",
        "cohorts": summaries,
        "cohort_count": len(cohorts),
        "pooled": pooled,
        "pair_coverage_passed": all(row["pair_coverage_passed"] for row in summaries),
        "missing_pairs": sum(row["missing_pairs"] for row in summaries),
        "extra_pairs": sum(row["extra_pairs"] for row in summaries),
        "abstained_pairs": sum(row["abstained_pairs"] for row in summaries),
        "calibration": _calibration(confidence, correctness)
        if complete_probability_coverage else None,
        "calibration_available": complete_probability_coverage,
        "bootstrap": bootstrap,
        "release_authorized": False,
        "scope": "Scoring only; caller must establish prospective freezing, independence, biological coverage and acceptance criteria.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--bootstrap-replicates", type=int, default=2000)
    parser.add_argument("--bootstrap-seed", type=int, default=20261120)
    args = parser.parse_args()
    original = args.input.read_bytes()
    data = json.loads(original)
    report = evaluate_cohorts(
        data["cohorts"], data["labels"],
        bootstrap_replicates=args.bootstrap_replicates,
        bootstrap_seed=args.bootstrap_seed,
    )
    if args.input.read_bytes() != original:
        raise RuntimeError("Scoring input changed during evaluation")
    report["input_sha256"] = hashlib.sha256(original).hexdigest()
    report["scorer_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("x") as handle:
        handle.write(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"cohorts": len(data["cohorts"]), "pairs": report["pooled"]["pairs"]}))


if __name__ == "__main__":
    main()
