from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--expected", type=int, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    records = [json.loads(path.read_text()) for path in sorted(args.input_dir.glob("task_*.json"))]
    if len(records) != args.expected:
        raise ValueError(f"expected {args.expected} records, found {len(records)}")
    methods = list(records[0]["estimates"])
    summaries = {}
    for method in methods:
        errors = np.asarray([record["absolute_error"][method] for record in records], dtype=float)
        summaries[method] = {
            "mean_absolute_error": float(np.mean(errors)),
            "median_absolute_error": float(np.median(errors)),
            "quantiles": {
                "0.5": float(np.quantile(errors, 0.5)),
                "0.9": float(np.quantile(errors, 0.9)),
                "0.95": float(np.quantile(errors, 0.95)),
                "1.0": float(np.max(errors)),
            },
            "within_10kb": float(np.mean(errors <= 10_000)),
            "within_25kb": float(np.mean(errors <= 25_000)),
            "within_50kb": float(np.mean(errors <= 50_000)),
        }
    selected = min(methods, key=lambda method: summaries[method]["median_absolute_error"])
    output = {
        "state": "COMPLETED",
        "replicates": len(records),
        "background_model": "Zigzag_1S14",
        "parameter_cell_design": {
            "selection_coefficients": sorted(
                {record["truth"]["selection_coefficient"] for record in records}
            ),
            "mutation_generations": sorted(
                {record["truth"]["mutation_generation"] for record in records}
            ),
            "minimum_end_frequencies": sorted(
                {record["truth"]["minimum_end_frequency"] for record in records}
            ),
        },
        "methods": summaries,
        "selected_by_median_absolute_error": selected,
        "records": records,
        "submission_count": 0,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()

