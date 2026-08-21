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
    keys = ("windows", "sfs", "targets", "ancestral_ne", "nuisance")
    missing = {
        path: [key for key in keys if key not in chunk]
        for path, chunk in zip(paths, chunks)
        if any(key not in chunk for key in keys)
    }
    if missing:
        raise ValueError(f"simulation chunks have incompatible schemas: {missing}")
    merged = {key: np.concatenate([chunk[key] for chunk in chunks]) for key in keys}
    merged["n_haplotypes"] = chunks[0]["n_haplotypes"]
    merged["nuisance_names"] = chunks[0]["nuisance_names"]
    merged["chunk_seeds"] = np.asarray([int(chunk["seed"]) for chunk in chunks])
    for chunk in chunks[1:]:
        if int(chunk["n_haplotypes"]) != int(merged["n_haplotypes"]):
            raise ValueError("simulation chunks disagree on n_haplotypes")
        if not np.array_equal(chunk["nuisance_names"], merged["nuisance_names"]):
            raise ValueError("simulation chunks disagree on nuisance schema")
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output, **merged)
    for chunk in chunks:
        chunk.close()


if __name__ == "__main__":
    main()
