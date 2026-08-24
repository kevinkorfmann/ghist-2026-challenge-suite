from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np


def genome_length(path: Path) -> float:
    bounds: dict[str, list[float]] = {}
    with path.open() as handle:
        for line in handle:
            chromosome, _, cm, _ = line.split()[:4]
            value = float(cm)
            if chromosome not in bounds:
                bounds[chromosome] = [value, value]
            else:
                bounds[chromosome][0] = min(bounds[chromosome][0], value)
                bounds[chromosome][1] = max(bounds[chromosome][1], value)
    return sum(end - start for start, end in bounds.values())


def read_segments(path: Path) -> dict[tuple[str, str], list[tuple[str, float]]]:
    pairs: dict[tuple[str, str], list[tuple[str, float]]] = defaultdict(list)
    with gzip.open(path, "rt") as handle:
        for line in handle:
            individual1, _, individual2, _, chromosome, _, _, length = line.split()
            key = tuple(sorted((individual1, individual2)))
            pairs[key].append((chromosome, float(length)))
    return pairs


def find(parent: dict[str, str], item: str) -> str:
    parent.setdefault(item, item)
    if parent[item] != item:
        parent[item] = find(parent, parent[item])
    return parent[item]


def union(parent: dict[str, str], left: str, right: str) -> None:
    root_left = find(parent, left)
    root_right = find(parent, right)
    if root_left != root_right:
        parent[root_right] = root_left


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ibd", type=Path, required=True)
    parser.add_argument("--hbd", type=Path, required=True)
    parser.add_argument("--map", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()

    genome_cm = genome_length(args.map)
    pairs = read_segments(args.ibd)
    rows = []
    parent: dict[str, str] = {}
    for (individual1, individual2), segments in pairs.items():
        lengths = np.asarray([length for _, length in segments])
        chromosomes = {chromosome for chromosome, _ in segments}
        record = {
            "indID1": individual1,
            "indID2": individual2,
            "segments": int(lengths.size),
            "chromosomes": len(chromosomes),
            "total_cm": float(np.sum(lengths)),
            "relatedness": float(np.sum(lengths)) / (2.0 * genome_cm),
            "mean_cm": float(np.mean(lengths)),
            "median_cm": float(np.median(lengths)),
            "max_cm": float(np.max(lengths)),
            "q90_cm": float(np.quantile(lengths, 0.90)),
        }
        for threshold in (5, 10, 20, 40):
            selected = lengths[lengths >= threshold]
            record[f"segments_ge_{threshold}"] = int(selected.size)
            record[f"total_cm_ge_{threshold}"] = float(np.sum(selected))
            record[f"relatedness_ge_{threshold}"] = float(np.sum(selected)) / (2.0 * genome_cm)
        long_relatedness = float(record["relatedness_ge_10"])
        record["degree_ge_10"] = (
            max(1, int(round(-math.log2(long_relatedness)))) if long_relatedness > 0 else 99
        )
        if long_relatedness >= 0.003 or float(record["max_cm"]) >= 20.0:
            union(parent, individual1, individual2)
        rows.append(record)

    hbd_totals: dict[str, float] = defaultdict(float)
    hbd_counts: dict[str, int] = defaultdict(int)
    with gzip.open(args.hbd, "rt") as handle:
        for line in handle:
            individual, _, _, _, _, _, _, length = line.split()
            hbd_totals[individual] += float(length)
            hbd_counts[individual] += 1

    args.out_dir.mkdir(parents=True, exist_ok=True)
    pair_path = args.out_dir / "pair_features.csv"
    with pair_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    hbd_path = args.out_dir / "hbd_features.csv"
    with hbd_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(("individual", "hbd_segments", "hbd_total_cm", "inbreeding_fraction"))
        for individual in sorted(hbd_totals):
            writer.writerow(
                (
                    individual,
                    hbd_counts[individual],
                    hbd_totals[individual],
                    hbd_totals[individual] / genome_cm,
                )
            )

    components: dict[str, list[str]] = defaultdict(list)
    for individual in parent:
        components[find(parent, individual)].append(individual)
    summary = {
        "state": "COMPLETED",
        "genome_cm": genome_cm,
        "pairs_with_ibd": len(rows),
        "relationship_graph_components": sorted(
            (sorted(component) for component in components.values()), key=len, reverse=True
        ),
        "component_sizes": sorted((len(component) for component in components.values()), reverse=True),
        "hbd_individuals": len(hbd_totals),
        "pair_features": str(pair_path),
        "hbd_features": str(hbd_path),
        "submission_count": 0,
    }
    (args.out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
