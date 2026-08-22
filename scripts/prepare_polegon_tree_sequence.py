from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import tskit


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    tree_sequence = tskit.load(args.input)
    tables = tree_sequence.dump_tables()
    before = np.asarray(tables.mutations.time).copy()
    tables.mutations.time = np.full(len(before), tskit.UNKNOWN_TIME)
    prepared = tables.tree_sequence()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    prepared.dump(args.output)
    report = {
        "input": str(args.input),
        "output": str(args.output),
        "samples": prepared.num_samples,
        "trees": prepared.num_trees,
        "sites": prepared.num_sites,
        "mutations": prepared.num_mutations,
        "explicit_mutation_times_removed": int(np.sum(~np.isnan(before))),
        "all_output_mutation_times_unknown": bool(
            np.all(np.isnan(prepared.tables.mutations.time))
        ),
    }
    args.output.with_suffix(args.output.suffix + ".json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
