"""Select a development candidate; never label old simulations held out."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def summarize(records):
    if len(records) != 48 or len({r["seed"] for r in records}) != 48:
        raise ValueError("Require all 48 distinct historical development simulations")
    names = records[0]["predictions"]
    if not names:
        raise ValueError("Require nonempty development candidates")
    configurations = {
        name: {
            key: prediction[key] for key in ("window_snps", "minimum_minor_frequency")
        }
        for name, prediction in names.items()
    }
    for record in records:
        if set(record["predictions"]) != set(names) or set(
            record["absolute_errors"]
        ) != set(names):
            raise ValueError("Development candidate sets differ")
        for key in ("source_record_sha256", "source_tree_sha256"):
            digest = record.get(key, "")
            if len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest):
                raise ValueError("Development records require source artifact hashes")
        origin, length = record["coordinate_origin"], record["sequence_length"]
        truth = record["truth"]["position"] - origin
        if (
            not np.isfinite([origin, length, truth]).all()
            or origin < 0
            or not 0 <= truth < length
        ):
            raise ValueError("Invalid development coordinate contract")
        for name, prediction in record["predictions"].items():
            if any(
                prediction[key] != value for key, value in configurations[name].items()
            ):
                raise ValueError("Development candidate configurations differ")
            position, stored = prediction["position"], record["absolute_errors"][name]
            if (
                not np.isfinite([position, stored]).all()
                or not 0 <= position < length
                or stored < 0
            ):
                raise ValueError(
                    "Development errors/positions must be finite and valid"
                )
            if stored != abs(position - truth):
                raise ValueError(
                    "Stored development errors differ from predictions and truth"
                )
    results = {}
    for name in names:
        errors = np.asarray([r["absolute_errors"][name] for r in records])
        results[name] = {
            "median_absolute_error": float(np.median(errors)),
            "q90_absolute_error": float(np.quantile(errors, 0.9)),
            "q95_absolute_error": float(np.quantile(errors, 0.95)),
            "within_75kb": int(sum(errors < 75_000)),
            "within_250kb": int(sum(errors < 250_000)),
            "informative": sum(r["predictions"][name]["informative"] for r in records),
            "median_maximum_ties": float(
                np.median([r["predictions"][name]["maximum_ties"] for r in records])
            ),
            "objective_median_plus_quarter_q90": float(
                np.median(errors) + 0.25 * np.quantile(errors, 0.9)
            ),
        }
    selected = min(
        results, key=lambda x: results[x]["objective_median_plus_quarter_q90"]
    )
    return {
        "status": "DEVELOPMENT ONLY: no held-out validation, calibration or promotion claim",
        "replicates": len(records),
        "source": "all48 historical Q10 standing-variation diagnostics now repurposed for development",
        "selection_objective": "median absolute error plus 0.25 times 90th-percentile error",
        "selected_candidate": selected,
        "development_metrics": results,
        "selected_configuration": {
            key: records[0]["predictions"][selected][key]
            for key in ("window_snps", "minimum_minor_frequency")
        },
        "selected_needs_fresh_confirmation": True,
        "records": records,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = summarize(
        [
            json.loads(path.read_text())
            for path in sorted(args.input_dir.glob("task_*.json"))
        ]
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "records"}, indent=2))


if __name__ == "__main__":
    main()
