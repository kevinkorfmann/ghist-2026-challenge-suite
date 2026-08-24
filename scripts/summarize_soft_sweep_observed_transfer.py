from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def metrics(records: list[dict], method: str, half_width: int) -> dict:
    errors = np.asarray(
        [record["absolute_error"][method] for record in records], dtype=float
    )
    coverage = float(np.mean(errors <= half_width))
    return {
        "replicates": len(records),
        "median_absolute_error": float(np.median(errors)),
        "q90_absolute_error": float(np.quantile(errors, 0.9)),
        "coverage": coverage,
        "mean_competition_score": coverage / (1.0 + half_width / 100_000.0),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--recovery-summary", type=Path, required=True)
    parser.add_argument("--observed-audit", type=Path, required=True)
    parser.add_argument("--signal-fraction", type=float, default=0.6)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    recovery = json.loads(args.recovery_summary.read_text())
    observed = json.loads(args.observed_audit.read_text())
    observed_score = float(
        observed["scales"]["25000"]["top_windows"][0]["soft_composite_score"]
    )
    threshold = args.signal_fraction * observed_score
    evaluation = [
        record for record in recovery["records"] if record["seed"] - 75_226_000 >= 24
    ]
    method = recovery["selected_method"]
    half_width = recovery["interval_width_selection"]["selected_half_width"]
    signal_matched = [
        record
        for record in evaluation
        if record["diagnostics"]["peak_soft_composite_score"] >= threshold
    ]
    high_frequency = [
        record
        for record in evaluation
        if record["truth"]["minimum_end_frequency"] >= 0.8
    ]
    signal_metrics = metrics(signal_matched, method, half_width)
    frequency_metrics = metrics(high_frequency, method, half_width)
    gate = {
        "at_least_eight_signal_matched_replicates": len(signal_matched) >= 8,
        "signal_matched_median_at_most_100kb": signal_metrics[
            "median_absolute_error"
        ]
        <= 100_000,
        "signal_matched_q90_at_most_250kb": signal_metrics["q90_absolute_error"]
        <= 250_000,
        "signal_matched_score_at_least_0.4": signal_metrics[
            "mean_competition_score"
        ]
        >= 0.4,
        "high_frequency_q90_at_most_250kb": frequency_metrics[
            "q90_absolute_error"
        ]
        <= 250_000,
    }
    gate["passed"] = all(gate.values())
    output = {
        "state": "COMPLETED",
        "selected_method": method,
        "selected_half_width": half_width,
        "observed_soft_composite_score": observed_score,
        "signal_fraction": args.signal_fraction,
        "signal_threshold": threshold,
        "signal_matched_untouched_evaluation": signal_metrics,
        "high_frequency_untouched_evaluation": frequency_metrics,
        "all_regime_untouched_evaluation": recovery["untouched_evaluation"],
        "robustness_gate": gate,
        "submission_count": 0,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
