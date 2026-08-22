from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--truth", type=Path, required=True)
    parser.add_argument("--singer", type=Path, required=True)
    parser.add_argument("--polegon", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    truth = np.load(args.truth)["mean"]
    methods = {
        "singer": np.load(args.singer)["mean"],
        "singer_polegon": np.load(args.polegon)["mean"],
    }
    result = {}
    for name, estimate in methods.items():
        result[name] = {
            "root_mean_squared_error": float(np.sqrt(np.mean((estimate - truth) ** 2))),
            "pair_class_rmse": np.sqrt(np.mean((estimate - truth) ** 2, axis=1)).tolist(),
            "pair_class_correlation": [
                float(np.corrcoef(estimate[row], truth[row])[0, 1]) for row in range(3)
            ],
        }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
