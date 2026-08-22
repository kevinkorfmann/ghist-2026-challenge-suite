from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def error_metrics(errors: np.ndarray) -> dict:
    return {
        "replicates": int(errors.size),
        "median_absolute_error": float(np.median(errors)),
        "quantiles": {
            "0.9": float(np.quantile(errors, 0.9)),
            "0.95": float(np.quantile(errors, 0.95)),
            "1.0": float(np.max(errors)),
        },
        "within_50kb": float(np.mean(errors <= 50_000)),
        "within_100kb": float(np.mean(errors <= 100_000)),
    }


def method_estimates(records: list[dict]) -> dict[str, np.ndarray]:
    gap = np.asarray([record["estimates"]["largest_gap_midpoint"] for record in records])
    composite = np.asarray([record["estimates"]["composite"] for record in records])
    linkage = np.asarray([record["estimates"]["flank_ld_100kb"] for record in records])
    return {
        "largest_gap_midpoint": gap,
        "median_all_three": np.median(np.vstack([gap, composite, linkage]), axis=0),
    }


def summarize(records: list[dict]) -> dict:
    truth = np.asarray([record["truth"]["position"] for record in records], dtype=float)
    estimates = method_estimates(records)
    scales = np.asarray([record["slim_scaling_factor"] for record in records], dtype=float)
    by_scale = {}
    for scale in sorted(set(scales)):
        keep = scales == scale
        by_scale[str(int(scale))] = {
            method: error_metrics(np.abs(values[keep] - truth[keep]))
            for method, values in estimates.items()
        }
    overall = {
        method: error_metrics(np.abs(values - truth)) for method, values in estimates.items()
    }
    selected = "median_all_three"
    selected_scale_metrics = {
        scale: methods[selected] for scale, methods in by_scale.items()
    }
    gate = {
        "median_at_most_25kb_each_scale": all(
            metrics["median_absolute_error"] <= 25_000
            for metrics in selected_scale_metrics.values()
        ),
        "q90_at_most_100kb_each_scale": all(
            metrics["quantiles"]["0.9"] <= 100_000
            for metrics in selected_scale_metrics.values()
        ),
    }
    gate["passed"] = all(gate.values())
    return {
        "state": "COMPLETED",
        "design": "eight matched chromosome-7 parameter cells at scale factors 1, 2, 5, and 10",
        "replicates": len(records),
        "selected_method_locked_before_stress_test": selected,
        "by_scale": by_scale,
        "overall": overall,
        "robustness_gate": gate,
        "submission_count": 0,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--expected", type=int, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    records = [json.loads(path.read_text()) for path in sorted(args.input_dir.glob("task_*.json"))]
    if len(records) != args.expected:
        raise ValueError(f"expected {args.expected} records, found {len(records)}")
    output = summarize(records)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
