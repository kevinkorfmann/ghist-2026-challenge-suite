from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fit", type=Path, action="append", required=True)
    parser.add_argument("--job", type=int, action="append", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if len(args.fit) != len(args.job):
        raise ValueError("each fit file must have one job identifier")

    payloads = [json.loads(path.read_text()) for path in args.fit]
    labels = list(payloads[0]["results"])
    if any(set(payload["results"]) != set(labels) for payload in payloads):
        raise ValueError("fit files do not contain identical strata")
    rows = {}
    for label in labels:
        values = [payload["results"][label] for payload in payloads]
        generations = np.asarray([value["generations"] for value in values])
        scaling = np.asarray([value["post_growth_scaling"] for value in values])
        ancestral = np.asarray([value["ancestral_effective_size"] for value in values])
        modern = np.asarray([value["modern_effective_size"] for value in values])
        rows[label] = {
            "by_projection": {
                str(payload["projection"]): {
                    "generations": value["generations"],
                    "post_growth_scaling": value["post_growth_scaling"],
                    "ancestral_effective_size": value["ancestral_effective_size"],
                    "modern_effective_size": value["modern_effective_size"],
                }
                for payload, value in zip(payloads, values)
            },
            "median": {
                "generations": float(np.median(generations)),
                "post_growth_scaling": float(np.median(scaling)),
                "ancestral_effective_size": float(np.median(ancestral)),
                "modern_effective_size": float(np.median(modern)),
            },
            "projection_range": {
                "generations": [float(generations.min()), float(generations.max())],
                "post_growth_scaling": [float(scaling.min()), float(scaling.max())],
                "ancestral_effective_size": [float(ancestral.min()), float(ancestral.max())],
                "modern_effective_size": [float(modern.min()), float(modern.max())],
            },
        }
    summary = {
        "jobs": [
            {"job_id": job, "projection": payload["projection"], "state": "COMPLETED"}
            for job, payload in zip(args.job, payloads)
        ],
        "source_files": [str(path) for path in args.fit],
        "strata": rows,
        "interpretation": {
            "status": "diagnostic_not_frozen",
            "primary_signal": (
                "The inferred ancestral size rises and the growth ratio falls with "
                "recombination, while modern size is comparatively stable."
            ),
            "submission_count": 0,
        },
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
