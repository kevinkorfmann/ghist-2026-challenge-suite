from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def describe(values: np.ndarray) -> dict:
    return {
        "mean": float(np.mean(values)),
        "standard_deviation": float(np.std(values, ddof=1)),
        "quantiles": {
            key: float(value)
            for key, value in zip(
                ("0.025", "0.25", "0.5", "0.75", "0.975"),
                np.quantile(values, (0.025, 0.25, 0.5, 0.75, 0.975)),
            )
        },
        "minimum": float(np.min(values)),
        "maximum": float(np.max(values)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--job", type=int, required=True)
    args = parser.parse_args()
    files = sorted(args.input_dir.glob("task_*.json"))
    if not files:
        raise ValueError("no bootstrap task files found")
    payloads = [json.loads(path.read_text()) for path in files]
    results = [row for payload in payloads for row in payload["results"]]
    replicate_ids = sorted(int(row["replicate"]) for row in results)
    if replicate_ids != list(range(len(results))):
        raise ValueError("bootstrap replicate identifiers are incomplete or duplicated")
    generations = np.asarray([row["generations"] for row in results])
    scaling = np.asarray([row["post_growth_scaling"] for row in results])
    summary = {
        "job_id": args.job,
        "state": "COMPLETED",
        "task_files": len(files),
        "replicates": len(results),
        "projection": payloads[0]["projection"],
        "block_size": payloads[0]["block_size"],
        "full_blocks_by_contig": payloads[0]["full_blocks_by_contig"],
        "terminal_partial_blocks_excluded_from_resampling": True,
        "generations": describe(generations),
        "post_growth_scaling": describe(scaling),
        "correlation": float(np.corrcoef(generations, scaling)[0, 1]),
        "provisional_hypothesis": {
            "generations": 5000.0,
            "post_growth_scaling": 2.0,
            "inside_marginal_95_percent_intervals": {
                "generations": bool(
                    np.quantile(generations, 0.025)
                    <= 5000.0
                    <= np.quantile(generations, 0.975)
                ),
                "post_growth_scaling": bool(
                    np.quantile(scaling, 0.025)
                    <= 2.0
                    <= np.quantile(scaling, 0.975)
                ),
            },
        },
        "submission_count": 0,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
