from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path

import numpy as np


PARAMETER_NAMES = (
    "mainland_population_size",
    "island_population_size",
    "generations_since_split",
    "generations_since_migration",
    "migration_rate",
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--observed", type=float, nargs=5, required=True)
    args = parser.parse_args()
    paths = [Path(path) for path in sorted(glob.glob(args.inputs))]
    if not paths:
        raise FileNotFoundError(f"no bootstrap results match {args.inputs!r}")
    records = [json.loads(path.read_text()) for path in paths]
    estimates = np.asarray(
        [
            (
                record["population1_effective_size"],
                record["population2_effective_size"],
                record["generations_since_split"],
                record["generations_since_contact"],
                record["migration_rate"],
            )
            for record in records
        ],
        dtype=np.float64,
    )
    truth = np.asarray([record["truth"] for record in records], dtype=np.float64)
    ratios = estimates / truth
    log_ratios = np.log(ratios)
    observed = np.asarray(args.observed, dtype=np.float64)
    median_ratio = np.exp(np.median(log_ratios, axis=0))
    ratio_interval = np.quantile(ratios, (0.025, 0.975), axis=0)
    corrected = observed / median_ratio
    corrected_interval = np.stack(
        (observed / ratio_interval[1], observed / ratio_interval[0]), axis=0
    )
    scores = np.sqrt(np.sum((ratios - 1.0) ** 2, axis=1))
    result = {
        "replicates": len(records),
        "parameter_names": PARAMETER_NAMES,
        "median_estimate_to_truth_ratio": median_ratio.tolist(),
        "mean_relative_error": np.mean(ratios - 1.0, axis=0).tolist(),
        "median_absolute_relative_error": np.median(
            np.abs(ratios - 1.0), axis=0
        ).tolist(),
        "ratio_interval_95": ratio_interval.tolist(),
        "competition_rrmse": {
            "mean": float(scores.mean()),
            "median": float(np.median(scores)),
            "interval_95": np.quantile(scores, (0.025, 0.975)).tolist(),
        },
        "observed_raw_fit": dict(zip(PARAMETER_NAMES, observed.tolist())),
        "observed_bias_corrected": dict(zip(PARAMETER_NAMES, corrected.tolist())),
        "observed_sampling_interval_95": {
            name: interval.tolist()
            for name, interval in zip(PARAMETER_NAMES, corrected_interval.T)
        },
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
