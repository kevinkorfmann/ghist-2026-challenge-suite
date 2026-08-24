from __future__ import annotations

import csv
import re
from collections import defaultdict
from pathlib import Path
from typing import TypeAlias

import numpy as np


LABEL_PATTERN = re.compile(r"^(.+?)(\d+)_g")
SexMap: TypeAlias = dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]


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
    return {
        chromosome: tuple(np.asarray(column, dtype=float) for column in columns)
        for chromosome, columns in values.items()
    }


def read_segment_rows(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        first = handle.readline()
        handle.seek(0)
        if first.startswith("indID1,"):
            return list(csv.DictReader(handle))
        rows = []
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) != 9:
                raise ValueError(f"Expected 9 fields in {path}, found {len(fields)}")
            rows.append(
                dict(
                    zip(
                        (
                            "indID1",
                            "indID2",
                            "chr",
                            "start_bp",
                            "end_bp",
                            "type",
                            "start_cM",
                            "end_cM",
                            "length",
                        ),
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


def pair_features(
    rows: list[dict[str, str]], genome_cm: float, sex_map: SexMap | None = None
) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[(row["indID1"], row["indID2"])].append(row)

    output: list[dict[str, object]] = []
    for (individual1, individual2), segments in grouped.items():
        lengths1 = [float(row["length"]) for row in segments if row["type"] == "IBD1"]
        lengths2 = [float(row["length"]) for row in segments if row["type"] == "IBD2"]
        male_lengths: list[float] = []
        female_lengths: list[float] = []
        if sex_map is not None:
            for row in segments:
                positions, male_cm, female_cm = sex_map[row["chr"]]
                bounds = np.asarray([float(row["start_bp"]), float(row["end_bp"])])
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
        if sex_map is not None:
            record.update(_summaries(male_lengths, "male_map"))
            record.update(_summaries(female_lengths, "female_map"))
            differences = [female - male for male, female in zip(male_lengths, female_lengths, strict=True)]
            record.update(_summaries(differences, "female_minus_male"))
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
    return values, feature_names
