from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def describe(values: np.ndarray) -> dict:
    return {
        "mean": float(np.mean(values)),
        "median": float(np.median(values)),
        "root_mean_square": float(np.sqrt(np.mean(values**2))),
        "quantiles": {
            key: float(value)
            for key, value in zip(
                ("0.025", "0.25", "0.75", "0.975"),
                np.quantile(values, (0.025, 0.25, 0.75, 0.975)),
            )
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--job", type=int, required=True)
    args = parser.parse_args()
    files = sorted(args.input_dir.glob("task_*.json"))
    rows = [json.loads(path.read_text()) for path in files]
    if len(rows) != 36:
        raise ValueError(f"expected 36 recovery fits, found {len(rows)}")
    generation_error = np.asarray(
        [row["relative_error"]["generations"] for row in rows]
    )
    scaling_error = np.asarray(
        [row["relative_error"]["post_growth_scaling"] for row in rows]
    )
    by_truth = {}
    for row in rows:
        truth = row["truth_in_analysis_units"]
        key = f"generations_{truth['generations']:g}|scaling_{truth['post_growth_scaling']:g}"
        by_truth.setdefault(key, []).append(row)
    grouped = {}
    for key, group in sorted(by_truth.items()):
        grouped[key] = {
            "replicates": len(group),
            "generation_relative_error": describe(
                np.asarray([row["relative_error"]["generations"] for row in group])
            ),
            "scaling_relative_error": describe(
                np.asarray(
                    [row["relative_error"]["post_growth_scaling"] for row in group]
                )
            ),
        }
    summary = {
        "job_id": args.job,
        "state": "COMPLETED",
        "replicates": len(rows),
        "generator": {
            "chromosome_arms": sorted(
                {row["generator"]["chromosome"] for row in rows}
            ),
            "stdpopsim": rows[0]["generator"]["stdpopsim"],
            "dfe": rows[0]["generator"]["dfe"],
            "genetic_map": rows[0]["generator"]["genetic_map"],
            "slim_scaling_factor": rows[0]["generator"]["slim_scaling_factor"],
            "slim_burn_in": rows[0]["generator"]["slim_burn_in"],
        },
        "overall": {
            "generation_relative_error": describe(generation_error),
            "scaling_relative_error": describe(scaling_error),
            "competition_rrmse": float(
                np.sqrt(np.mean(np.concatenate((generation_error, scaling_error)) ** 2))
            ),
        },
        "by_truth": grouped,
        "submission_count": 0,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
