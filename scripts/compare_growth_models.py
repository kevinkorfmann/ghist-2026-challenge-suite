from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from ghist_growth_bgs.selection import compare_growth_models


def main():
    parser = argparse.ArgumentParser(
        description="Compare growth histories on withheld chromosome arms"
    )
    parser.add_argument("--blocks", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--projection", type=int, default=40)
    parser.add_argument("--restarts", type=int, default=12)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--maxiter", type=int, default=500)
    parser.add_argument("--seed", type=int, default=20261120)
    args = parser.parse_args()
    with np.load(args.blocks) as data:
        result = compare_growth_models(
            data["spectra"],
            data["exposure_bp"],
            data["contig"],
            projection=args.projection,
            restarts=args.restarts,
            workers=args.workers,
            maxiter=args.maxiter,
            seed=args.seed,
        )
    result["source"] = str(args.blocks)
    result["source_sha256"] = hashlib.sha256(args.blocks.read_bytes()).hexdigest()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps(
            {
                key: value
                for key, value in result.items()
                if key not in {"folds", "full_data_fits"}
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
