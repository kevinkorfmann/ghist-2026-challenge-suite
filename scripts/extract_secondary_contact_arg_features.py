from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path

import numpy as np
import tskit

from ghist_secondary_contact.features import population_pair_hazards


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trees", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    paths = [Path(path) for path in sorted(glob.glob(args.trees))]
    if not paths:
        raise FileNotFoundError(f"no tree sequences match {args.trees!r}")
    hazards = []
    summaries = []
    for path in paths:
        tree_sequence = tskit.load(path)
        samples = tree_sequence.samples()
        if len(samples) != 134:
            raise ValueError(f"expected 134 haplotypes, found {len(samples)} in {path}")
        value = population_pair_hazards(
            tree_sequence,
            samples[:86],
            samples[86:],
        )
        hazards.append(value)
        summaries.append(
            {
                "path": str(path),
                "sequence_length": tree_sequence.sequence_length,
                "trees": tree_sequence.num_trees,
                "nodes": tree_sequence.num_nodes,
                "edges": tree_sequence.num_edges,
            }
        )
    array = np.asarray(hazards)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.out,
        hazards=array,
        mean=array.mean(axis=0),
        standard_deviation=array.std(axis=0),
    )
    report = {
        "draws": len(paths),
        "shape": list(array.shape),
        "tree_sequences": summaries,
    }
    args.out.with_suffix(args.out.suffix + ".json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
