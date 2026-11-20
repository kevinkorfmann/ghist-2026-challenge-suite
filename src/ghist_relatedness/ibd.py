from __future__ import annotations

import csv
import re
from collections import defaultdict
from pathlib import Path
from typing import TypeAlias

import numpy as np


LABEL_PATTERN = re.compile(r"^(.+?)(\d+)_g")
SexMap: TypeAlias = dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]
LENGTH_BINS = np.asarray(
    [0, 2, 4, 6, 8, 10, 15, 20, 30, 40, 60, 80, 120, np.inf], dtype=float
)


def read_sex_map(path: Path) -> SexMap:
    values: dict[str, list[list[float]]] = defaultdict(lambda: [[], [], []])
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chromosome, position, male_cm, female_cm = line.split()[:4]
            values[chromosome][0].append(float(position))
            values[chromosome][1].append(float(male_cm))
            values[chromosome][2].append(float(female_cm))
    result = {
        chromosome: tuple(np.asarray(column, dtype=float) for column in columns)
        for chromosome, columns in values.items()
    }
    if not result:
        raise ValueError("Empty genetic map")
    for chromosome, columns in result.items():
        if (
            not all(np.all(np.isfinite(column)) for column in columns)
            or np.any(np.diff(columns[0]) <= 0)
            or any(np.any(np.diff(column) < 0) for column in columns[1:])
        ):
            raise ValueError(f"Nonfinite or nonmonotonic genetic map: {chromosome}")
    return result


def read_segment_rows(path: Path) -> list[dict[str, str]]:
    columns = (
        "indID1", "indID2", "chr", "start_bp", "end_bp", "type",
        "start_cM", "end_cM", "length",
    )
    with path.open(encoding="utf-8-sig") as handle:
        lines = (line for line in handle if line.strip() and not line.startswith("#"))
        first = next(lines, "")
        if not first:
            raise ValueError("Empty IBD segment input")
        from itertools import chain

        lines = chain((first,), lines)
        if "," in first:
            reader = csv.DictReader(lines)
            if len(reader.fieldnames or []) != len(columns) or set(reader.fieldnames) != set(columns):
                raise ValueError(f"Expected unique IBD CSV columns: {columns}")
            rows = list(reader)
            if any(None in row or any(value is None or not value.strip() for value in row.values()) for row in rows):
                raise ValueError("Missing or extra IBD CSV fields")
            return rows
        rows = []
        for line in lines:
            fields = line.split()
            if len(fields) != 9:
                raise ValueError(f"Expected 9 fields in {path}, found {len(fields)}")
            rows.append(
                dict(
                    zip(
                        columns,
                        fields,
                        strict=True,
                    )
                )
            )
        return rows


def relationship_from_sample(sample: str) -> str:
    match = LABEL_PATTERN.match(sample)
    if match is None:
        raise ValueError(f"Cannot recover relationship label from {sample}")
    return match.group(1)


def _summaries(values: list[float], prefix: str) -> dict[str, float]:
    if not values:
        return {f"{prefix}_{name}": 0.0 for name in ("mean", "sd", "q25", "median", "q75", "q90", "max")}
    array = np.asarray(values, dtype=float)
    return {
        f"{prefix}_mean": float(np.mean(array)),
        f"{prefix}_sd": float(np.std(array)),
        f"{prefix}_q25": float(np.quantile(array, 0.25)),
        f"{prefix}_median": float(np.median(array)),
        f"{prefix}_q75": float(np.quantile(array, 0.75)),
        f"{prefix}_q90": float(np.quantile(array, 0.90)),
        f"{prefix}_max": float(np.max(array)),
    }


def _distribution_features(
    values: list[float], prefix: str, ordered: int = 20
) -> dict[str, float]:
    array = np.asarray(values, dtype=float)
    counts, _ = np.histogram(array, bins=LENGTH_BINS)
    output = {
        f"{prefix}_bin_{index}": float(count)
        for index, count in enumerate(counts)
    }
    descending = np.sort(array)[::-1]
    output.update(
        {
            f"{prefix}_ordered_{index}": float(descending[index])
            if index < descending.size
            else 0.0
            for index in range(ordered)
        }
    )
    return output


def pair_features(
    rows: list[dict[str, str]], genome_cm: float, sex_map: SexMap | None = None
) -> list[dict[str, object]]:
    if not np.isfinite(genome_cm) or genome_cm <= 0:
        raise ValueError("genome_cm must be positive and finite")
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        pair = tuple(sorted((row["indID1"], row["indID2"])))
        if not all(pair) or pair[0] == pair[1]:
            raise ValueError("IBD segments must identify two distinct individuals")
        if row["type"] not in {"IBD1", "IBD2"}:
            raise ValueError(f"Unsupported IBD segment type: {row['type']}")
        length = float(row["length"])
        if not np.isfinite(length) or length < 0:
            raise ValueError("IBD segment length must be finite and nonnegative")
        grouped[pair].append(row)

    output: list[dict[str, object]] = []
    for (individual1, individual2), segments in grouped.items():
        lengths1 = [float(row["length"]) for row in segments if row["type"] == "IBD1"]
        lengths2 = [float(row["length"]) for row in segments if row["type"] == "IBD2"]
        male_lengths: list[float] = []
        female_lengths: list[float] = []
        if sex_map is not None:
            for row in segments:
                if row["chr"] not in sex_map:
                    raise ValueError(f"Chromosome absent from genetic map: {row['chr']}")
                positions, male_cm, female_cm = sex_map[row["chr"]]
                bounds = np.asarray([float(row["start_bp"]), float(row["end_bp"])])
                if not np.all(np.isfinite(bounds)) or bounds[0] < 0 or bounds[1] < bounds[0]:
                    raise ValueError("Segment physical coordinates must be finite and ordered")
                male_bounds = np.interp(bounds, positions, male_cm)
                female_bounds = np.interp(bounds, positions, female_cm)
                male_lengths.append(float(male_bounds[1] - male_bounds[0]))
                female_lengths.append(float(female_bounds[1] - female_bounds[0]))
        total1 = float(sum(lengths1))
        total2 = float(sum(lengths2))
        chromosome_totals: dict[str, float] = defaultdict(float)
        chromosome_types: dict[str, set[str]] = defaultdict(set)
        boundary_starts = 0
        for row in segments:
            chromosome_totals[row["chr"]] += float(row["length"])
            chromosome_types[row["chr"]].add(row["type"])
            boundary_starts += float(row["start_cM"]) <= 1e-6
        chromosome_values = list(chromosome_totals.values())
        record: dict[str, object] = {
            "indID1": individual1,
            "indID2": individual2,
            "total_ibd1_cm": total1,
            "total_ibd2_cm": total2,
            "relatedness": (0.5 * total1 + total2) / genome_cm,
            "ibd2_share": total2 / max(total1 + total2, 1e-9),
            "n_ibd1": float(len(lengths1)),
            "n_ibd2": float(len(lengths2)),
            "n_segments": float(len(segments)),
            "n_chromosomes": float(len(chromosome_totals)),
            "n_ibd2_chromosomes": float(sum("IBD2" in types for types in chromosome_types.values())),
            "boundary_start_fraction": boundary_starts / len(segments),
            "mean_cm_per_chromosome": float(np.mean(chromosome_values)),
            "sd_cm_per_chromosome": float(np.std(chromosome_values)),
            "max_cm_per_chromosome": float(np.max(chromosome_values)),
        }
        record.update(_summaries(lengths1, "ibd1"))
        record.update(_summaries(lengths2, "ibd2"))
        record.update(_distribution_features(lengths1, "ibd1"))
        record.update(_distribution_features(lengths2, "ibd2"))
        record.update(_distribution_features(chromosome_values, "chromosome_total", ordered=22))
        if sex_map is not None:
            record.update(_summaries(male_lengths, "male_map"))
            record.update(_summaries(female_lengths, "female_map"))
            differences = [female - male for male, female in zip(male_lengths, female_lengths, strict=True)]
            record.update(_summaries(differences, "female_minus_male"))
            record.update(_distribution_features(male_lengths, "male_map"))
            record.update(_distribution_features(female_lengths, "female_map"))
            record["total_male_map_cm"] = float(sum(male_lengths))
            record["total_female_map_cm"] = float(sum(female_lengths))
            record["female_male_total_ratio"] = float(sum(female_lengths)) / max(
                float(sum(male_lengths)), 1e-9
            )
        output.append(record)
    return output


def matrix(records: list[dict[str, object]], feature_names: list[str] | None = None) -> tuple[np.ndarray, list[str]]:
    if not records:
        raise ValueError("No pair records")
    if feature_names is None:
        feature_names = sorted(name for name in records[0] if name not in {"indID1", "indID2", "label"})
    values = np.asarray([[float(record[name]) for name in feature_names] for record in records])
    if not np.all(np.isfinite(values)):
        raise ValueError("Nonfinite relatedness features")
    return values, feature_names
