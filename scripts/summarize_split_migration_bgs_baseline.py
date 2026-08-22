from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


PARAMETERS = (
    "AFR_population_size",
    "EUR_population_size",
    "generations_since_split",
    "migration_rate",
    "ancestral_population_size",
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("files", nargs="+", type=Path)
    parser.add_argument("--job-id", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    rows = [json.loads(path.read_text()) for path in args.files]
    fits = {"x".join(map(str, row["fit"]["projection"])): row["fit"] for row in rows}
    payload = {
        "state": "COMPLETED",
        "stratum": "noncoding_1kb_100kb",
        "physical_exposure_bp": rows[0]["physical_exposure_bp"],
        "job_id": args.job_id,
        "by_projection": fits,
        "projection_median": {
            key: float(np.median([row["fit"][key] for row in rows]))
            for key in PARAMETERS
        },
        "projection_range": {
            key: [
                float(min(row["fit"][key] for row in rows)),
                float(max(row["fit"][key] for row in rows)),
            ]
            for key in PARAMETERS
        },
        "metadata": rows[0]["metadata"],
        "submission_count": 0,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
