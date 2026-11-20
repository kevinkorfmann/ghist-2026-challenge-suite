from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from ghist_single_sweep.recovery import validate_blind_records


ENSEMBLES = {
    "median_soft_haplotype": ("soft_composite", "h12", "h12_h2"),
    "median_soft_linkage": ("soft_composite", "h12", "flank_ld_100kb"),
    "median_soft_four": ("soft_composite", "h12", "h12_h2", "flank_ld_100kb"),
}


def metrics(records: list[dict], method: str, half_width: int) -> dict:
    if not records:
        return {
            "replicates": 0,
            "median_absolute_error": None,
            "q90_absolute_error": None,
            "coverage": None,
            "mean_competition_score": None,
        }
    names = ENSEMBLES.get(method, (method,))
    centers = np.asarray(
        [np.median([r["estimates"][name] for name in names]) for r in records]
    )
    truth = np.asarray([r["truth"]["position"] for r in records])
    errors = np.abs(centers - truth)
    coverage = float(
        np.mean((centers - half_width <= truth) & (truth < centers + half_width))
    )
    return {
        "replicates": len(records),
        "median_absolute_error": float(np.median(errors)),
        "q90_absolute_error": float(np.quantile(errors, 0.9)),
        "coverage": coverage,
        "mean_competition_score": coverage / (1 + half_width / 100_000),
    }


def summarize_transfer(
    recovery: dict, observed: dict, signal_fraction=0.6, window_size=100_000
):
    validate_blind_records(recovery.get("records", []))
    if not 0 < signal_fraction <= 1:
        raise ValueError("Signal fraction must lie in (0,1]")
    if any(r.get("scan_window_size") != window_size for r in recovery["records"]):
        raise ValueError(
            "Observed and simulated score window sizes must match explicitly; regenerate missing scale metadata"
        )
    try:
        observed_score = float(
            observed["scales"][str(window_size)]["top_windows"][0][
                "soft_composite_score"
            ]
        )
    except (KeyError, IndexError) as error:
        raise ValueError(
            "Observed audit lacks the matching scan window size"
        ) from error
    if not np.isfinite(observed_score) or observed_score <= 0:
        raise ValueError("Observed soft sweep score is not informative")
    threshold = signal_fraction * observed_score
    evaluation = [r for r in recovery["records"] if r["seed"] - 75_226_000 >= 24]
    method = recovery["selected_method"]
    half_width = recovery["interval_width_selection"]["selected_half_width"]
    if half_width <= 0:
        raise ValueError("Interval half-width must be positive")
    signal_matched = [
        r
        for r in evaluation
        if r["diagnostics"]["peak_soft_composite_score"] >= threshold
    ]
    high_frequency = [
        r for r in evaluation if r["truth"]["minimum_end_frequency"] >= 0.8
    ]
    signal_metrics, frequency_metrics = (
        metrics(rows, method, half_width) for rows in (signal_matched, high_frequency)
    )
    gate = {
        "at_least_eight_signal_matched_replicates": len(signal_matched) >= 8,
        "signal_matched_median_at_most_100kb": bool(signal_matched)
        and signal_metrics["median_absolute_error"] <= 100_000,
        "signal_matched_q90_at_most_250kb": bool(signal_matched)
        and signal_metrics["q90_absolute_error"] <= 250_000,
        "signal_matched_score_at_least_0.4": bool(signal_matched)
        and signal_metrics["mean_competition_score"] >= 0.4,
        "high_frequency_q90_at_most_250kb": bool(high_frequency)
        and frequency_metrics["q90_absolute_error"] <= 250_000,
        "not_repurposed_or_retrospective": not any(
            r.get("reanalysis_of_legacy_simulation") or r.get("usage") == "development"
            for r in recovery["records"]
        ),
    }
    gate["passed"] = all(gate.values())
    return {
        "state": "COMPLETED",
        "evaluation_status": recovery.get(
            "evaluation_status", "provenance requires independent audit"
        ),
        "selected_method": method,
        "selected_half_width": half_width,
        "score_window_size": window_size,
        "observed_soft_composite_score": observed_score,
        "signal_fraction": signal_fraction,
        "signal_threshold": threshold,
        "signal_matched_evaluation": signal_metrics,
        "high_frequency_evaluation": frequency_metrics,
        "all_regime_evaluation": recovery["untouched_evaluation"],
        "robustness_gate": gate,
        "submission_count": 0,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--recovery-summary", type=Path, required=True)
    parser.add_argument("--observed-audit", type=Path, required=True)
    parser.add_argument("--signal-fraction", type=float, default=0.6)
    parser.add_argument("--window-size", type=int, default=100_000)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    output = summarize_transfer(
        json.loads(args.recovery_summary.read_text()),
        json.loads(args.observed_audit.read_text()),
        args.signal_fraction,
        args.window_size,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
