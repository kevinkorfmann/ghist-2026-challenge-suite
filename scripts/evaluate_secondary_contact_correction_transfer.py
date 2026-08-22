from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path

import numpy as np

from ghist_secondary_contact.calibration import (
    LogCorrection,
    apply_log_correction,
    competition_scores,
)


ESTIMATE_KEYS = (
    "population1_effective_size",
    "population2_effective_size",
    "generations_since_split",
    "generations_since_contact",
    "migration_rate",
)


def labelled_pattern(value: str) -> tuple[str, str]:
    label, separator, pattern = value.partition("=")
    if not separator or not label or not pattern:
        raise argparse.ArgumentTypeError("expected LABEL=GLOB")
    return label, pattern


def summary(predictions: np.ndarray, truth: np.ndarray) -> dict:
    scores = competition_scores(predictions, truth)
    relative = predictions / truth - 1.0
    return {
        "mean_competition_rrmse": float(scores.mean()),
        "median_competition_rrmse": float(np.median(scores)),
        "competition_rrmse_interval_95": np.quantile(scores, (0.025, 0.975)).tolist(),
        "median_absolute_relative_error": dict(
            zip(ESTIMATE_KEYS, np.median(np.abs(relative), axis=0).tolist())
        ),
        "mean_relative_error": dict(zip(ESTIMATE_KEYS, relative.mean(axis=0).tolist())),
        "replicate_scores": scores.tolist(),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--correction", type=Path, required=True)
    parser.add_argument("--condition", action="append", type=labelled_pattern, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    correction_source = json.loads(args.correction.read_text())
    serialized = correction_source["final_correction"]
    if serialized is None or correction_source["accepted_prediction_method"] == "identity":
        correction = None
    else:
        correction = LogCorrection(
            center=np.asarray(serialized["center"], dtype=np.float64),
            scale=np.asarray(serialized["scale"], dtype=np.float64),
            coefficients=np.asarray(serialized["coefficients"], dtype=np.float64),
        )

    conditions = {}
    for label, pattern in args.condition:
        paths = sorted(Path(path) for path in glob.glob(pattern))
        if not paths:
            raise FileNotFoundError(f"no results match {pattern!r}")
        records = [json.loads(path.read_text()) for path in paths]
        estimates = np.asarray(
            [[record[key] for key in ESTIMATE_KEYS] for record in records],
            dtype=np.float64,
        )
        truth = np.asarray([record["truth"] for record in records], dtype=np.float64)
        corrected = apply_log_correction(correction, estimates)
        raw_scores = competition_scores(estimates, truth)
        corrected_scores = competition_scores(corrected, truth)
        difference = corrected_scores - raw_scores
        conditions[label] = {
            "pattern": pattern,
            "replicates": len(records),
            "raw": summary(estimates, truth),
            "corrected": summary(corrected, truth),
            "corrected_minus_raw": {
                "mean": float(difference.mean()),
                "median": float(np.median(difference)),
                "interval_95": np.quantile(difference, (0.025, 0.975)).tolist(),
                "corrected_wins": int(np.sum(difference < 0)),
                "ties": int(np.sum(difference == 0)),
                "raw_wins": int(np.sum(difference > 0)),
                "replicate_differences": difference.tolist(),
            },
        }
    result = {
        "correction": str(args.correction),
        "accepted_prediction_method": correction_source["accepted_prediction_method"],
        "conditions": conditions,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
