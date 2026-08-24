from __future__ import annotations

import argparse
import copy
import json
from collections import defaultdict
from pathlib import Path

import joblib
import numpy as np

from fit_human_family_relationships import (
    BIN_EDGES,
    MEIOSES,
    ORDERED,
    background_correct,
    matrix,
    relationship_label,
    summarize_raw_lengths,
)


def ordered_lengths(record: dict[str, object]) -> list[float]:
    return [
        float(record[f"ordered_length_{rank + 1}"])
        for rank in range(ORDERED)
        if float(record[f"ordered_length_{rank + 1}"]) > 0
    ]


def read_pedsim_with_lengths(
    path: Path, genome_cm: float
) -> tuple[list[dict[str, object]], dict[tuple[str, str], list[float]]]:
    grouped: dict[tuple[str, str], list[tuple[str, float]]] = defaultdict(list)
    labels: dict[tuple[str, str], str] = {}
    with path.open() as handle:
        for line in handle:
            left, right, chromosome, _, _, state, _, _, length = line.split()
            pair = tuple(sorted((left, right)))
            labels[pair] = relationship_label(left)
            grouped[pair].append((chromosome, float(length)))
            if state == "IBD2":
                grouped[pair].append((chromosome, float(length)))
    records = []
    lengths = {}
    for pair in sorted(grouped):
        record: dict[str, object] = {
            "indID1": pair[0],
            "indID2": pair[1],
            "label": labels[pair],
        }
        record.update(summarize_raw_lengths(grouped[pair], genome_cm))
        records.append(record)
        lengths[pair] = [value for _, value in grouped[pair]]
    return records, lengths


def add_background(
    records: list[dict[str, object]],
    source_lengths: dict[tuple[str, str], list[float]],
    rng: np.random.Generator,
    genome_cm: float,
) -> list[dict[str, object]]:
    augmented = []
    individual_effects: dict[str, float] = {}
    samples = {
        str(record[key])
        for record in records
        for key in ("indID1", "indID2")
    }
    for sample in samples:
        individual_effects[sample] = float(rng.lognormal(mean=0.0, sigma=0.25))
    for record in records:
        effect = np.sqrt(
            individual_effects[str(record["indID1"])]
            * individual_effects[str(record["indID2"])]
        )
        count = int(rng.poisson(3.0 * effect))
        background = 20.0 + np.minimum(rng.exponential(10.0, size=count), 60.0)
        pair = tuple(sorted((str(record["indID1"]), str(record["indID2"]))))
        lengths = source_lengths[pair] + background.tolist()
        output = {
            "indID1": record["indID1"],
            "indID2": record["indID2"],
            "label": record["label"],
        }
        output.update(summarize_raw_lengths([("0", value) for value in lengths], genome_cm))
        augmented.append(output)
    return augmented


def trim_shortest_background_shape(
    records: list[dict[str, object]], count_offset: float, minimum_cm: int
) -> None:
    remove_count = max(0, round(count_offset))
    for record in records:
        lengths = ordered_lengths(record)
        eligible = sorted((value for value in lengths if value >= minimum_cm))
        removed = eligible[:remove_count]
        retained = list(lengths)
        for value in removed:
            retained.remove(value)
        retained.sort(reverse=True)
        record["max_cm"] = retained[0] if retained else 0.0
        removed_counts, _ = np.histogram(
            np.asarray(removed), bins=np.asarray(BIN_EDGES)
        )
        for index, count in enumerate(removed_counts):
            key = f"length_bin_{index}"
            record[key] = max(0.0, float(record[key]) - float(count))
        for rank in range(ORDERED):
            record[f"ordered_length_{rank + 1}"] = (
                retained[rank] if rank < len(retained) else 0.0
            )


def accuracy(model, records: list[dict[str, object]], names: list[str]) -> tuple[float, float]:
    truth = [str(record["label"]) for record in records]
    predicted = model.predict(matrix(records, names))
    exact = float(np.mean([left == right for left, right in zip(truth, predicted, strict=True)]))
    degree = float(
        np.mean(
            [
                MEIOSES[left] == MEIOSES[str(right)]
                for left, right in zip(truth, predicted, strict=True)
            ]
        )
    )
    return exact, degree


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--evaluation", type=Path, required=True)
    parser.add_argument("--replicates", type=int, default=12)
    parser.add_argument("--genome-cm", type=float, default=3346.298892538)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    payload = joblib.load(args.model)
    model, names = payload["model"], payload["feature_names"]
    minimum_cm = min(
        int(name.rsplit("_", 1)[1])
        for name in names
        if name.startswith("segments_ge_")
    )
    reference, source_lengths = read_pedsim_with_lengths(
        args.evaluation, args.genome_cm
    )
    baseline = accuracy(model, reference, names)
    rows = []
    for replicate in range(args.replicates):
        augmented = add_background(
            reference,
            source_lengths,
            np.random.default_rng(121226000 + replicate),
            args.genome_cm,
        )
        for quantile in (0.25, 0.5, 0.75):
            corrected = copy.deepcopy(augmented)
            offsets = background_correct(
                corrected, reference, quantile, args.genome_cm
            )
            for trim in (False, True):
                evaluated = copy.deepcopy(corrected)
                if trim:
                    trim_shortest_background_shape(
                        evaluated,
                        float(offsets[f"segments_ge_{minimum_cm}"]),
                        minimum_cm,
                    )
                exact, degree = accuracy(model, evaluated, names)
                rows.append(
                    {
                        "replicate": replicate,
                        "quantile": quantile,
                        "trim_shortest_shape": trim,
                        "exact_accuracy": exact,
                        "meioses_accuracy": degree,
                    }
                )
    aggregates = []
    for quantile in (0.25, 0.5, 0.75):
        for trim in (False, True):
            selected = [
                row
                for row in rows
                if row["quantile"] == quantile and row["trim_shortest_shape"] == trim
            ]
            aggregates.append(
                {
                    "quantile": quantile,
                    "trim_shortest_shape": trim,
                    "mean_exact_accuracy": float(
                        np.mean([row["exact_accuracy"] for row in selected])
                    ),
                    "mean_meioses_accuracy": float(
                        np.mean([row["meioses_accuracy"] for row in selected])
                    ),
                }
            )
    summary = {
        "state": "COMPLETED",
        "design": "held-out human pairs with heterogeneous added long-segment background",
        "replicates": args.replicates,
        "baseline_exact_accuracy": baseline[0],
        "baseline_meioses_accuracy": baseline[1],
        "aggregates": aggregates,
        "runs": rows,
        "submission_count": 0,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({key: value for key, value in summary.items() if key != "runs"}, indent=2))


if __name__ == "__main__":
    main()
