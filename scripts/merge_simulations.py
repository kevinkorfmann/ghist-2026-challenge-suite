from __future__ import annotations

import argparse
import glob
from pathlib import Path

import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    paths = sorted(glob.glob(args.input))
    if not paths:
        raise FileNotFoundError(args.input)
    chunks = [np.load(path) for path in paths]
    keys = ("windows", "sfs", "targets", "ancestral_ne")
    merged = {key: np.concatenate([chunk[key] for chunk in chunks]) for key in keys}
    merged["n_haplotypes"] = chunks[0]["n_haplotypes"]
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output, **merged)


if __name__ == "__main__":
    main()
