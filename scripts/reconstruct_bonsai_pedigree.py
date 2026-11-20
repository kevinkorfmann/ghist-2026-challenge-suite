from __future__ import annotations

import argparse
import gzip
import json
from collections import defaultdict
from itertools import combinations
from pathlib import Path
from typing import NamedTuple

import numpy as np


class HaplotypeSegment(NamedTuple):
    start: int
    end: int
    haplotype1: int
    haplotype2: int


def read_map(path: Path) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    values: dict[str, list[list[float]]] = defaultdict(lambda: [[], []])
    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.split()
            if not fields:
                continue
            if len(fields) != 4:
                raise ValueError(
                    f"{path}:{line_number}: expected four PLINK map fields"
                )
            chromosome, _, cm, bp = fields
            try:
                position, genetic_position = float(bp), float(cm)
            except ValueError as error:
                raise ValueError(
                    f"{path}:{line_number}: invalid map coordinate"
                ) from error
            if not np.isfinite(position) or position <= 0 or not position.is_integer():
                raise ValueError(
                    f"{path}:{line_number}: map bp must be a positive integer"
                )
            if not np.isfinite(genetic_position) or genetic_position < 0:
                raise ValueError(
                    f"{path}:{line_number}: map cM must be finite and nonnegative"
                )
            values[chromosome][0].append(position)
            values[chromosome][1].append(genetic_position)
    if not values:
        raise ValueError(f"{path}: empty genetic map")
    result = {
        chromosome: (np.asarray(columns[0]), np.asarray(columns[1]))
        for chromosome, columns in values.items()
    }
    for chromosome, (positions, cm) in result.items():
        if len(positions) < 2 or np.any(np.diff(positions) <= 0):
            raise ValueError(
                f"{path}: chromosome {chromosome} needs increasing, unique map bp"
            )
        if np.any(np.diff(cm) < 0):
            raise ValueError(f"{path}: chromosome {chromosome} has decreasing map cM")
    return result


def numeric_id(sample: str) -> int:
    try:
        return int(sample.rsplit("_", 1)[1]) + 1
    except (IndexError, ValueError) as error:
        raise ValueError(
            f"Expected an identifier ending in an integer: {sample}"
        ) from error


def read_haplotype_segments(
    path: Path,
    allowed_samples: set[str] | None = None,
) -> tuple[dict[tuple[str, str, str], list[HaplotypeSegment]], set[str]]:
    grouped: dict[tuple[str, str, str], list[HaplotypeSegment]] = defaultdict(list)
    samples: set[str] = set()
    with gzip.open(path, "rt") as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.split()
            if not fields:
                continue
            if len(fields) != 8:
                raise ValueError(f"{path}:{line_number}: expected eight hap-IBD fields")
            left, hap1, right, hap2, chromosome, start, end, length = fields
            try:
                segment = HaplotypeSegment(int(start), int(end), int(hap1), int(hap2))
                _validate_segment(segment)
                reported_cm = float(length)
                if not np.isfinite(reported_cm) or reported_cm < 0:
                    raise ValueError("reported cM must be finite and nonnegative")
                if left == right:
                    raise ValueError("within-sample HBD is not between-sample IBD")
            except ValueError as error:
                raise ValueError(f"{path}:{line_number}: {error}") from error
            if allowed_samples is not None and (
                left not in allowed_samples or right not in allowed_samples
            ):
                continue
            if left > right:
                left, right = right, left
                segment = HaplotypeSegment(
                    segment.start, segment.end, segment.haplotype2, segment.haplotype1
                )
            grouped[(left, right, chromosome)].append(segment)
            samples.update((left, right))
    return grouped, samples


def _validate_segment(segment: tuple[int, int, int, int]) -> None:
    if len(segment) != 4:
        raise ValueError("IBD intervals require (start, end, haplotype1, haplotype2)")
    if any(
        isinstance(value, bool) or not isinstance(value, (int, np.integer))
        for value in segment
    ):
        raise ValueError("IBD coordinates and haplotype identifiers must be integers")
    start, end, hap1, hap2 = segment
    if start <= 0 or end <= start:
        raise ValueError("IBD coordinates require 0 < start < end")
    if hap1 not in (1, 2) or hap2 not in (1, 2):
        raise ValueError("hap-IBD haplotype identifiers must be 1 or 2")


def coverage_intervals(
    intervals: list[tuple[int, int, int, int]],
) -> list[tuple[int, int, bool]]:
    """Union haplotype edges into IBD1/IBD2 spans for one ordered sample pair.

    Full sharing requires a two-edge matching: 1->1 with 2->2, or 1->2 with
    2->1. Repeated edges and edges sharing a haplotype cannot create IBD2.
    Counts retain an edge until its last overlapping record ends. Physical
    endpoints are unchanged; cM lengths integrate map(end) - map(start), so
    endpoints meeting at one marker do not create a positive-length overlap.
    Ambiguous legacy (start, end) records are deliberately rejected.
    """
    events: dict[int, dict[tuple[int, int], int]] = defaultdict(
        lambda: defaultdict(int)
    )
    for segment in intervals:
        _validate_segment(segment)
        start, end, hap1, hap2 = segment
        events[start][(hap1, hap2)] += 1
        events[end][(hap1, hap2)] -= 1
    output: list[tuple[int, int, bool]] = []
    active: dict[tuple[int, int], int] = defaultdict(int)
    previous: int | None = None
    for position in sorted(events):
        if previous is not None and position > previous and any(active.values()):
            is_full = bool(
                (active[(1, 1)] and active[(2, 2)])
                or (active[(1, 2)] and active[(2, 1)])
            )
            if output and output[-1][1] == previous and output[-1][2] == is_full:
                output[-1] = (output[-1][0], position, is_full)
            else:
                output.append((previous, position, is_full))
        for edge, change in events[position].items():
            active[edge] += change
        previous = position
    return output


def convert_segments(
    path: Path,
    map_path: Path,
    threshold_cm: float,
    allowed_samples: set[str] | None = None,
    cm_scale: float = 1.0,
) -> tuple[list[list[object]], list[dict[str, object]], dict[int, str]]:
    if not np.isfinite(threshold_cm) or threshold_cm < 0:
        raise ValueError("threshold_cm must be finite and nonnegative")
    if not np.isfinite(cm_scale) or cm_scale <= 0:
        raise ValueError("cm_scale must be finite and positive")
    genetic_map = read_map(map_path)
    grouped, samples = read_haplotype_segments(path, allowed_samples)
    if allowed_samples is not None:
        samples.update(allowed_samples)
    id_to_sample = {numeric_id(sample): sample for sample in samples}
    if len(id_to_sample) != len(samples):
        raise ValueError("Sample identifiers are not uniquely mapped")
    converted: list[list[object]] = []
    for (sample1, sample2, chromosome), intervals in grouped.items():
        if chromosome not in genetic_map:
            raise ValueError(
                f"IBD chromosome {chromosome} is missing from the genetic map"
            )
        positions, cm = genetic_map[chromosome]
        for start, end, is_full in coverage_intervals(intervals):
            if start < positions[0] or end > positions[-1]:
                raise ValueError(
                    f"IBD interval {chromosome}:{start}-{end} exceeds genetic map bounds"
                )
            bounds = np.interp(np.asarray([start, end], dtype=float), positions, cm)
            length_cm = float(bounds[1] - bounds[0])
            # Thresholds use the original species' map; scale only retained lengths.
            if length_cm > 0 and length_cm >= threshold_cm:
                scaled_length = length_cm * cm_scale
                if not np.isfinite(scaled_length):
                    raise ValueError("scaled segment cM must be finite")
                converted.append(
                    [
                        numeric_id(sample1),
                        numeric_id(sample2),
                        chromosome,
                        start,
                        end,
                        is_full,
                        scaled_length,
                    ]
                )
    bio_info = [
        {"genotype_id": identifier, "age": None, "sex": None}
        for identifier in sorted(id_to_sample)
    ]
    return converted, bio_info, id_to_sample


def serializable_topology(
    up_dict: dict[object, list[object]],
) -> dict[str, list[object]]:
    return {str(identifier): values for identifier, values in up_dict.items()}


def segment_summaries(
    segments: list[list[object]], identifiers: list[int]
) -> dict[frozenset[int], dict[str, float]]:
    summaries = {
        frozenset((left, right)): {
            "total_half": 0.0,
            "total_full": 0.0,
            "num_half": 0,
            "max_seg_cm": 0.0,
        }
        for left, right in combinations(identifiers, 2)
    }
    for identifier1, identifier2, _, _, _, is_full, length_cm in segments:
        summary = summaries[frozenset((int(identifier1), int(identifier2)))]
        key = "total_full" if is_full else "total_half"
        summary[key] += float(length_cm)
        summary["num_half"] += int(not is_full)
        summary["max_seg_cm"] = max(summary["max_seg_cm"], float(length_cm))
    return summaries


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ibd", type=Path, required=True)
    parser.add_argument("--map", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--threshold-cm", type=float, default=10.0)
    parser.add_argument("--max-add-degree", type=int, default=4)
    parser.add_argument("--max-rel-types", type=int, default=3)
    parser.add_argument("--cm-scale", type=float, default=1.0)
    args = parser.parse_args()

    from bonsaitree import bonsai

    segments, bio_info, id_to_sample = convert_segments(
        args.ibd, args.map, args.threshold_cm, cm_scale=args.cm_scale
    )

    # The challenge provides no ages or sexes. The implementation supports missing
    # metadata internally, so only relax its public input guard; likelihoods then
    # integrate over the missing values instead of inventing observations.
    bonsai.validate_bio = lambda bio, idx: None
    result = bonsai.build_pedigree(
        segments,
        bio_info,
        ibd_stats=segment_summaries(segments, sorted(id_to_sample)),
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
    print(
        json.dumps(
            {
                key: value
                for key, value in summary.items()
                if key not in {"topology", "relationships"}
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
