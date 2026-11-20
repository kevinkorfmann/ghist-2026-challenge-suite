from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from ghist_single_sweep.recovery import validate_blind_records


HALF_WIDTHS = np.asarray(
    [25_000, 35_000, 50_000, 60_000, 75_000, 100_000, 125_000, 150_000, 200_000]
)


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


def summarize(records: list[dict], seed_origin: int = 75_226_000) -> dict:
    validate_blind_records(records)
    if len(records) != 48:
        raise ValueError(f"expected 48 records, found {len(records)}")
    truth = np.asarray([record["truth"]["position"] for record in records], dtype=float)
    raw = {
        name: np.asarray([record["estimates"][name] for record in records], dtype=float)
        for name in records[0]["estimates"]
    }
    methods = {
        **raw,
        "median_soft_haplotype": np.median(
            np.vstack([raw["soft_composite"], raw["h12"], raw["h12_h2"]]), axis=0
        ),
        "median_soft_linkage": np.median(
            np.vstack([raw["soft_composite"], raw["h12"], raw["flank_ld_100kb"]]),
            axis=0,
        ),
        "median_soft_four": np.median(
            np.vstack(
                [
                    raw["soft_composite"],
                    raw["h12"],
                    raw["h12_h2"],
                    raw["flank_ld_100kb"],
                ]
            ),
            axis=0,
        ),
    }
    indices = np.asarray([record["seed"] - seed_origin for record in records], dtype=int)
    calibration = indices < 24
    evaluation = ~calibration
    selection_metrics = {
        name: error_metrics(np.abs(values[calibration] - truth[calibration]))
        for name, values in methods.items()
    }
    selection_objective = {
        name: metrics["median_absolute_error"] + 0.25 * metrics["quantiles"]["0.9"]
        for name, metrics in selection_metrics.items()
    }
    selected = min(selection_objective, key=selection_objective.get)
    selected_errors = np.abs(methods[selected] - truth)
    width_scores = {
        str(int(width)): float(
            np.mean(selected_errors[calibration] <= width) / (1.0 + width / 100_000.0)
        )
        for width in HALF_WIDTHS
    }
    selected_width = int(max(width_scores, key=width_scores.get))
    evaluation_metrics = error_metrics(selected_errors[evaluation])
    evaluation_coverage = float(np.mean(selected_errors[evaluation] <= selected_width))
    evaluation_score = evaluation_coverage / (1.0 + selected_width / 100_000.0)
    gate = {
        "evaluation_median_at_most_100kb": evaluation_metrics["median_absolute_error"] <= 100_000,
        "evaluation_q90_at_most_250kb": evaluation_metrics["quantiles"]["0.9"] <= 250_000,
        "evaluation_mean_competition_score_at_least_0.4": evaluation_score >= 0.4,
    }
    gate["passed"] = all(gate.values())
    return {
        "state": "COMPLETED",
        "evaluation_status": "retrospective reanalysis of previously evaluated simulations" if any(record.get("reanalysis_of_legacy_simulation") for record in records) else "seed-heldout evaluation; provenance requires independent audit",
        "design": "48 chromosome-map standing-variation sweeps; first 24 method calibration and last 24 seed-heldout evaluation",
        "replicates": len(records),
        "calibration_replicates": int(calibration.sum()),
        "evaluation_replicates": int(evaluation.sum()),
        "selection_metrics": selection_metrics,
        "selection_objective_median_plus_quarter_q90": selection_objective,
        "selected_method": selected,
        "interval_width_selection": {
            "candidate_half_widths": HALF_WIDTHS.tolist(),
            "calibration_scores": width_scores,
            "selected_half_width": selected_width,
        },
        "untouched_evaluation": {
            "metrics": evaluation_metrics,
            "candidate_coverage": evaluation_coverage,
            "mean_competition_score": evaluation_score,
        },
        "robustness_gate": gate,
        "records": records,
        "submission_count": 0,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--expected", type=int, default=48)
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
