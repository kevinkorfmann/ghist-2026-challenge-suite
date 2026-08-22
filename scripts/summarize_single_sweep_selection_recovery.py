from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def error_metrics(errors: np.ndarray) -> dict:
    return {
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--raisd-dir", type=Path)
    parser.add_argument("--expected", type=int, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    records = [json.loads(path.read_text()) for path in sorted(args.input_dir.glob("task_*.json"))]
    if len(records) != args.expected:
        raise ValueError(f"expected {args.expected} records, found {len(records)}")
    candidate_estimates: dict[str, np.ndarray] = {
        method: np.asarray([record["estimates"][method] for record in records], dtype=float)
        for method in records[0]["estimates"]
    }
    gap = candidate_estimates["largest_gap_midpoint"]
    composite = candidate_estimates["composite"]
    linkage = candidate_estimates["flank_ld_100kb"]
    reference = None
    if args.raisd_dir is not None:
        reference_records = [
            json.loads((args.raisd_dir / f"task_{task_index}.json").read_text())
            for task_index in range(args.expected)
        ]
        reference = np.asarray(
            [record["maximum"]["position"] for record in reference_records], dtype=float
        )
        candidate_estimates["raisd"] = reference
    candidate_estimates.update(
        {
            "median_all_three": np.median(np.vstack([gap, composite, linkage]), axis=0),
            "midpoint_gap_composite": (gap + composite) / 2.0,
            "midpoint_gap_linkage": (gap + linkage) / 2.0,
            "midpoint_composite_linkage": (composite + linkage) / 2.0,
            "agreement_gated_gap_composite": np.where(
                np.abs(gap - composite) <= 100_000,
                (gap + composite) / 2.0,
                gap,
            ),
        }
    )
    if reference is not None:
        candidate_estimates.update(
            {
                "median_with_raisd": np.median(
                    np.vstack([gap, composite, linkage, reference]), axis=0
                ),
                "median_gap_composite_raisd": np.median(
                    np.vstack([gap, composite, reference]), axis=0
                ),
            }
        )
    truth = np.asarray([record["truth"]["position"] for record in records], dtype=float)
    task_indices = np.asarray([record["seed"] - 75_222_000 for record in records])
    partition_code = (task_indices * 17 + 3) % 5
    selection_mask = np.isin(partition_code, [0, 1])
    evaluation_mask = ~selection_mask
    summaries = {
        method: error_metrics(np.abs(estimates - truth))
        for method, estimates in candidate_estimates.items()
    }
    selection_metrics = {
        method: error_metrics(np.abs(estimates[selection_mask] - truth[selection_mask]))
        for method, estimates in candidate_estimates.items()
    }
    objective = {
        method: metrics["median_absolute_error"] + 0.25 * metrics["quantiles"]["0.9"]
        for method, metrics in selection_metrics.items()
    }
    selected = min(objective, key=objective.get)
    evaluation_metrics = {
        method: error_metrics(np.abs(estimates[evaluation_mask] - truth[evaluation_mask]))
        for method, estimates in candidate_estimates.items()
    }
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
        "held_out_method_gate": {
            "partition": "parameter cells assigned by a fixed seed hash",
            "selection_cells": int(selection_mask.sum()),
            "evaluation_cells": int(evaluation_mask.sum()),
            "selection_metrics": selection_metrics,
            "selection_objective_median_plus_quarter_q90": objective,
            "selected": selected,
            "evaluation_metrics": evaluation_metrics,
        },
        "records": records,
        "submission_count": 0,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
