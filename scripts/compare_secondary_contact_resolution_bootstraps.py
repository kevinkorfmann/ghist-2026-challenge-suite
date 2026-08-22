from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path

import numpy as np


PARAMETER_KEYS = (
    "population1_effective_size",
    "population2_effective_size",
    "generations_since_split",
    "generations_since_contact",
    "migration_rate",
)
SCALED_PARAMETER_KEYS = (
    "population1_size_ratio",
    "population2_size_ratio",
    "split_duration_2N",
    "contact_fraction",
    "migration_2N",
)
LOWER_BOUNDS = np.asarray((0.01, 0.01, 1e-5, 0.002, 1e-6))
UPPER_BOUNDS = np.asarray((100.0, 100.0, 10.0, 0.98, 100.0))


def labelled_pattern(value: str) -> tuple[str, str]:
    label, separator, pattern = value.partition("=")
    if not separator or not label or not pattern:
        raise argparse.ArgumentTypeError("expected LABEL=GLOB")
    return label, pattern


def load_resolution(pattern: str) -> dict[int, dict]:
    records: dict[int, dict] = {}
    for filename in sorted(glob.glob(pattern)):
        record = json.loads(Path(filename).read_text())
        replicate = int(record["replicate"])
        if replicate in records:
            raise ValueError(f"duplicate replicate {replicate} for {pattern!r}")
        records[replicate] = record
    if not records:
        raise FileNotFoundError(f"no results match {pattern!r}")
    return records


def percentile_interval(values: np.ndarray) -> list[float]:
    return np.quantile(values, (0.025, 0.975)).tolist()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--resolution", action="append", type=labelled_pattern, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--bootstrap-draws", type=int, default=20_000)
    parser.add_argument("--seed", type=int, default=20260822)
    args = parser.parse_args()
    if len(args.resolution) < 2:
        parser.error("at least two --resolution values are required")

    loaded = {label: load_resolution(pattern) for label, pattern in args.resolution}
    replicate_sets = {label: set(records) for label, records in loaded.items()}
    common = set.intersection(*replicate_sets.values())
    if any(indices != common for indices in replicate_sets.values()):
        details = {label: sorted(indices) for label, indices in replicate_sets.items()}
        raise ValueError(f"replicate sets differ: {details}")
    replicate_ids = sorted(common)

    truth_by_replicate: dict[int, np.ndarray] = {}
    for replicate in replicate_ids:
        truths = [
            np.asarray(loaded[label][replicate]["truth"], dtype=np.float64)
            for label in loaded
        ]
        if any(not np.array_equal(truths[0], truth) for truth in truths[1:]):
            raise ValueError(f"truth differs across resolutions for replicate {replicate}")
        truth_by_replicate[replicate] = truths[0]

    scores: dict[str, np.ndarray] = {}
    relative_errors: dict[str, np.ndarray] = {}
    summaries = {}
    for label, records in loaded.items():
        estimates = np.asarray(
            [[records[index][key] for key in PARAMETER_KEYS] for index in replicate_ids],
            dtype=np.float64,
        )
        truth = np.asarray([truth_by_replicate[index] for index in replicate_ids])
        errors = estimates / truth - 1.0
        score = np.sqrt(np.sum(errors**2, axis=1))
        relative_errors[label] = errors
        scores[label] = score
        boundary_replicates = []
        for index in replicate_ids:
            scaled = np.asarray(
                [records[index]["scaled_parameters"][key] for key in SCALED_PARAMETER_KEYS]
            )
            near_boundary = (scaled <= LOWER_BOUNDS * 1.001) | (
                scaled >= UPPER_BOUNDS / 1.001
            )
            if np.any(near_boundary):
                boundary_replicates.append(
                    {
                        "replicate": index,
                        "parameters": [
                            key
                            for key, flagged in zip(SCALED_PARAMETER_KEYS, near_boundary)
                            if flagged
                        ],
                    }
                )
        summaries[label] = {
            "mean_competition_rrmse": float(score.mean()),
            "median_competition_rrmse": float(np.median(score)),
            "competition_rrmse_interval_95": percentile_interval(score),
            "median_absolute_relative_error": dict(
                zip(PARAMETER_KEYS, np.median(np.abs(errors), axis=0).tolist())
            ),
            "replicate_scores": dict(zip(map(str, replicate_ids), score.tolist())),
            "boundary_replicates": boundary_replicates,
        }

    rng = np.random.default_rng(args.seed)
    draw_indices = rng.integers(0, len(replicate_ids), (args.bootstrap_draws, len(replicate_ids)))
    pairwise = {}
    labels = list(loaded)
    for left_index, left in enumerate(labels):
        for right in labels[left_index + 1 :]:
            delta = scores[left] - scores[right]
            mean_draws = delta[draw_indices].mean(axis=1)
            pairwise[f"{left}_minus_{right}"] = {
                "mean_score_difference": float(delta.mean()),
                "median_score_difference": float(np.median(delta)),
                "mean_difference_bootstrap_interval_95": percentile_interval(mean_draws),
                "left_wins": int(np.sum(delta < 0.0)),
                "ties": int(np.sum(delta == 0.0)),
                "right_wins": int(np.sum(delta > 0.0)),
                "replicate_differences": dict(zip(map(str, replicate_ids), delta.tolist())),
            }

    ranking = sorted(
        labels,
        key=lambda label: (
            summaries[label]["mean_competition_rrmse"],
            summaries[label]["median_competition_rrmse"],
        ),
    )
    selection_valid = not any(
        summary["boundary_replicates"] for summary in summaries.values()
    )
    result = {
        "replicates": replicate_ids,
        "parameter_keys": PARAMETER_KEYS,
        "selection_rule": "lowest mean competition RRMSE, then lowest median competition RRMSE",
        "ranking": ranking,
        "selection_valid": selection_valid,
        "selected_resolution": ranking[0] if selection_valid else None,
        "summaries": summaries,
        "paired_comparisons": pairwise,
        "bootstrap_draws": args.bootstrap_draws,
        "seed": args.seed,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
