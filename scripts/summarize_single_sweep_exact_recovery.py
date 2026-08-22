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


def summarize(records: list[dict]) -> dict:
    truth = np.asarray([record["truth"]["position"] for record in records], dtype=float)
    raw = np.asarray(
        [
            [
                record["estimates"]["largest_gap_midpoint"],
                record["estimates"]["composite"],
                record["estimates"]["flank_ld_100kb"],
            ]
            for record in records
        ],
        dtype=float,
    )
    estimate = np.median(raw, axis=1)
    errors = np.abs(estimate - truth)
    spread = np.ptp(raw, axis=1)
    consensus = spread <= 50_000
    seeds = np.asarray([record["seed"] for record in records], dtype=int)
    calibration = seeds < 75_225_064
    evaluation = ~calibration
    candidate_widths = np.asarray(
        [25_000, 35_000, 40_000, 45_000, 50_000, 55_000, 60_000, 75_000, 100_000]
    )
    calibration_consensus = calibration & consensus
    width_scores = {
        str(int(width)): float(
            np.mean((errors[calibration_consensus] <= width) / (1.0 + width / 100_000.0))
        )
        for width in candidate_widths
    }
    half_width = int(max(width_scores, key=width_scores.get))
    evaluation_consensus = evaluation & consensus
    overall = error_metrics(errors)
    evaluation_overall = error_metrics(errors[evaluation])
    evaluation_consensus_metrics = (
        error_metrics(errors[evaluation_consensus]) if evaluation_consensus.any() else None
    )
    hit_score = 1.0 / (1.0 + half_width / 100_000.0)
    gate = {
        "calibration_has_at_least_eight_consensus_replicates": int(calibration_consensus.sum()) >= 8,
        "evaluation_has_at_least_six_consensus_replicates": int(evaluation_consensus.sum()) >= 6,
        "evaluation_overall_median_at_most_75kb": evaluation_overall["median_absolute_error"] <= 75_000,
        "evaluation_overall_q90_at_most_250kb": evaluation_overall["quantiles"]["0.9"] <= 250_000,
        "evaluation_consensus_coverage_at_least_75pct": bool(
            evaluation_consensus.any()
            and np.mean(errors[evaluation_consensus] <= half_width) >= 0.75
        ),
        "evaluation_consensus_q90_within_selected_width": bool(
            evaluation_consensus.any()
            and np.quantile(errors[evaluation_consensus], 0.9) <= half_width
        ),
    }
    gate["passed"] = all(gate.values())
    return {
        "state": "COMPLETED",
        "design": f"{len(records)} untouched exact-scale chromosome-7 map simulations across eight parameter cells",
        "replicates": len(records),
        "selected_method_locked_before_campaign": "median_all_three",
        "observed_estimator_spread": 46_000,
        "consensus_definition": "three primary estimators span at most 50 kb",
        "calibration_replicates": int(calibration.sum()),
        "evaluation_replicates": int(evaluation.sum()),
        "calibration_consensus_replicates": int(calibration_consensus.sum()),
        "evaluation_consensus_replicates": int(evaluation_consensus.sum()),
        "overall": overall,
        "interval_width_selection": {
            "candidate_half_widths": candidate_widths.tolist(),
            "calibration_consensus_scores": width_scores,
            "selected_half_width": half_width,
        },
        "untouched_evaluation": {
            "overall": evaluation_overall,
            "observed_like_consensus": evaluation_consensus_metrics,
            "consensus_candidate_coverage": float(
                np.mean(errors[evaluation_consensus] <= half_width)
            )
            if evaluation_consensus.any()
            else None,
        },
        "candidate_half_width": half_width,
        "candidate_hit_score": hit_score,
        "robustness_gate": gate,
        "records": records,
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
