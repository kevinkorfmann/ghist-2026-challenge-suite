from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from scripts.summarize_single_sweep_scale_stress import error_metrics


def summarize(records: list[dict], half_width: int = 35_000) -> dict:
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
    overall = error_metrics(errors)
    consensus_metrics = error_metrics(errors[consensus]) if consensus.any() else None
    hit_score = 1.0 / (1.0 + half_width / 100_000.0)
    gate = {
        "overall_median_at_most_75kb": overall["median_absolute_error"] <= 75_000,
        "overall_q90_at_most_250kb": overall["quantiles"]["0.9"] <= 250_000,
        "at_least_eight_observed_like_consensus_replicates": int(consensus.sum()) >= 8,
        "consensus_35kb_coverage_at_least_80pct": bool(
            consensus.any() and np.mean(errors[consensus] <= half_width) >= 0.8
        ),
        "consensus_q90_at_most_50kb": bool(
            consensus.any() and np.quantile(errors[consensus], 0.9) <= 50_000
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
        "consensus_replicates": int(consensus.sum()),
        "overall": overall,
        "consensus": consensus_metrics,
        "candidate_half_width": half_width,
        "candidate_hit_score": hit_score,
        "consensus_candidate_coverage": float(np.mean(errors[consensus] <= half_width))
        if consensus.any()
        else None,
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
