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
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--task-dir", type=Path, required=True)
    parser.add_argument("--block-summary", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    files = sorted(
        args.task_dir.glob("task_*.json"),
        key=lambda path: int(path.stem.rsplit("_", 1)[1]),
    )
    if len(files) != 64:
        raise ValueError(f"expected 64 bootstrap tasks, found {len(files)}")
    rows = [json.loads(path.read_text()) for path in files]
    values = {key: np.asarray([row["fit"][key] for row in rows]) for key in PARAMETERS}
    summary = {
        key: {
            "mean": float(np.mean(column)),
            "standard_deviation": float(np.std(column, ddof=1)),
            "quantiles": dict(
                zip(
                    ("0.025", "0.25", "0.5", "0.75", "0.975"),
                    np.quantile(column, (0.025, 0.25, 0.5, 0.75, 0.975)).tolist(),
                )
            ),
        }
        for key, column in values.items()
    }
    payload = {
        "state": "COMPLETED",
        "replicates": len(rows),
        "block_design": json.loads(args.block_summary.read_text()),
        "parameters": summary,
        "submission_count": 0,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()

