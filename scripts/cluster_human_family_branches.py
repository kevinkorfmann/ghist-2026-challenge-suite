from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


def find(parent: dict[str, str], sample: str) -> str:
    if parent[sample] != sample:
        parent[sample] = find(parent, parent[sample])
    return parent[sample]


def union(parent: dict[str, str], left: str, right: str) -> None:
    left_root, right_root = find(parent, left), find(parent, right)
    if left_root != right_root:
        parent[right_root] = left_root


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pairs", type=Path, required=True)
    parser.add_argument("--threshold", type=float, default=0.12)
    parser.add_argument("--field", default="relatedness_ge_20")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    with args.pairs.open() as handle:
        rows = list(csv.DictReader(handle))
    samples = sorted(
        {sample for row in rows for sample in (row["indID1"], row["indID2"])}
    )
    parent = {sample: sample for sample in samples}
    retained_edges = 0
    for row in rows:
        if float(row[args.field]) >= args.threshold:
            union(parent, row["indID1"], row["indID2"])
            retained_edges += 1
    grouped: dict[str, list[str]] = defaultdict(list)
    for sample in samples:
        grouped[find(parent, sample)].append(sample)
    families = sorted((sorted(group) for group in grouped.values()), key=lambda x: (-len(x), x))
    summary = {
        "state": "COMPLETED",
        "field": args.field,
        "threshold": args.threshold,
        "retained_edges": retained_edges,
        "families": families,
        "family_sizes": [len(family) for family in families],
        "samples": len(samples),
        "submission_count": 0,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({key: value for key, value in summary.items() if key != "families"}, indent=2))


if __name__ == "__main__":
    main()
