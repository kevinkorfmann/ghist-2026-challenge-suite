from __future__ import annotations

import argparse
import gzip
import json
from collections import defaultdict
from pathlib import Path

import numpy as np


def read_map(path: Path) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    values: dict[str, list[list[float]]] = defaultdict(lambda: [[], []])
    with path.open() as handle:
        for line in handle:
            chromosome, _, cm, bp = line.split()[:4]
            values[chromosome][0].append(float(bp))
            values[chromosome][1].append(float(cm))
    return {
        chromosome: (np.asarray(columns[0]), np.asarray(columns[1]))
        for chromosome, columns in values.items()
    }


def numeric_id(sample: str) -> int:
    try:
        return int(sample.rsplit("_", 1)[1]) + 1
    except (IndexError, ValueError) as error:
        raise ValueError(f"Expected an identifier ending in an integer: {sample}") from error


def read_haplotype_segments(
    path: Path,
) -> tuple[dict[tuple[str, str, str], list[tuple[int, int]]], set[str]]:
    grouped: dict[tuple[str, str, str], list[tuple[int, int]]] = defaultdict(list)
    samples: set[str] = set()
    with gzip.open(path, "rt") as handle:
        for line in handle:
            left, _, right, _, chromosome, start, end, _ = line.split()
            individual1, individual2 = sorted((left, right))
            grouped[(individual1, individual2, chromosome)].append((int(start), int(end)))
            samples.update((left, right))
    return grouped, samples


def coverage_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int, bool]]:
    events: dict[int, int] = defaultdict(int)
    for start, end in intervals:
        events[start] += 1
        events[end] -= 1
    output: list[tuple[int, int, bool]] = []
    coverage = 0
    previous: int | None = None
    for position in sorted(events):
        if previous is not None and position > previous and coverage > 0:
            is_full = coverage >= 2
            if output and output[-1][1] == previous and output[-1][2] == is_full:
                output[-1] = (output[-1][0], position, is_full)
            else:
                output.append((previous, position, is_full))
        coverage += events[position]
        previous = position
    return output


def convert_segments(
    path: Path,
    map_path: Path,
    threshold_cm: float,
) -> tuple[list[list[object]], list[dict[str, object]], dict[int, str]]:
    genetic_map = read_map(map_path)
    grouped, samples = read_haplotype_segments(path)
    id_to_sample = {numeric_id(sample): sample for sample in samples}
    if len(id_to_sample) != len(samples):
        raise ValueError("Sample identifiers are not uniquely mapped")
    converted: list[list[object]] = []
    for (sample1, sample2, chromosome), intervals in grouped.items():
        positions, cm = genetic_map[chromosome]
        for start, end, is_full in coverage_intervals(intervals):
            bounds = np.interp(np.asarray([start, end], dtype=float), positions, cm)
            length_cm = float(bounds[1] - bounds[0])
            if length_cm >= threshold_cm:
                converted.append(
                    [
                        numeric_id(sample1),
                        numeric_id(sample2),
                        chromosome,
                        start,
                        end,
                        is_full,
                        length_cm,
                    ]
                )
    bio_info = [
        {"genotype_id": identifier, "age": None, "sex": None}
        for identifier in sorted(id_to_sample)
    ]
    return converted, bio_info, id_to_sample


def serializable_topology(up_dict: dict[object, list[object]]) -> dict[str, list[object]]:
    return {str(identifier): values for identifier, values in up_dict.items()}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ibd", type=Path, required=True)
    parser.add_argument("--map", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--threshold-cm", type=float, default=10.0)
    parser.add_argument("--max-add-degree", type=int, default=4)
    parser.add_argument("--max-rel-types", type=int, default=3)
    args = parser.parse_args()

    from bonsaitree import bonsai

    segments, bio_info, id_to_sample = convert_segments(
        args.ibd, args.map, args.threshold_cm
    )

    # The challenge provides no ages or sexes. The implementation supports missing
    # metadata internally, so only relax its public input guard; likelihoods then
    # integrate over the missing values instead of inventing observations.
    bonsai.validate_bio = lambda bio, idx: None
    result = bonsai.build_pedigree(
        segments,
        bio_info,
        max_add_degree=args.max_add_degree,
        max_rel_append_types=args.max_rel_types,
        ibd_threshold=0,
        remove_distant_threshold=float("inf"),
    )
    pedigree = result["ped_obj"]
    genotyped = sorted(id_to_sample)
    relationship_counts: dict[str, int] = defaultdict(int)
    relationship_rows = []
    for offset, identifier1 in enumerate(genotyped):
        for identifier2 in genotyped[offset + 1 :]:
            relationship = pedigree.rels.get(identifier1, {}).get(identifier2)
            if relationship is None:
                continue
            relationship = tuple(relationship)
            relationship_counts[str(relationship)] += 1
            relationship_rows.append(
                {
                    "indID1": id_to_sample[identifier1],
                    "indID2": id_to_sample[identifier2],
                    "relationship_tuple": relationship,
                }
            )

    summary = {
        "state": "COMPLETED",
        "input_segments": len(segments),
        "individuals": len(genotyped),
        "placed_pairs": len(relationship_rows),
        "relationship_counts": dict(sorted(relationship_counts.items())),
        "topology": serializable_topology(pedigree.up_pedigree_dict),
        "relationships": relationship_rows,
        "metadata_policy": "age and sex unknown",
        "submission_count": 0,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({key: value for key, value in summary.items() if key not in {"topology", "relationships"}}, indent=2))


if __name__ == "__main__":
    main()
