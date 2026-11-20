from __future__ import annotations

import csv
import math
from itertools import combinations
from pathlib import Path

import yaml

from ghist_relatedness.labels import DOG_LABELS, HUMAN_LABELS, IBD_LABEL_MEIOSES

DEMOGRAPHY_FIELDS = {
    1: ("generations", "post_decline_fraction"),
    2: (
        "mainland_population_size",
        "island_population_size",
        "generations_since_split",
        "generations_since_migration",
        "migration_rate",
    ),
    3: ("generations", "post_growth_scaling"),
    4: (
        "AFR_population_size",
        "EUR_population_size",
        "generations_since_split",
        "migration_rate",
    ),
}


class UniqueKeyLoader(yaml.SafeLoader):
    """Reject duplicate YAML keys instead of silently retaining the last value."""


def _mapping(loader, node, deep=False):
    result = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=deep)
        if not isinstance(key, str) or key in result:
            raise ValueError(f"Duplicate or non-string YAML key: {key!r}")
        result[key] = loader.construct_object(value_node, deep=deep)
    return result


UniqueKeyLoader.add_constructor(
    yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _mapping
)


def validate_demography(path: Path, challenge: int) -> dict:
    data = yaml.load(path.read_text(), Loader=UniqueKeyLoader)
    fields = set(DEMOGRAPHY_FIELDS[challenge])
    if not isinstance(data, dict) or set(data) != fields:
        raise ValueError(f"Expected exactly these YAML fields: {sorted(fields)}")
    for name, value in data.items():
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
        ):
            raise ValueError(f"{name} must be a finite number")
        if value < 0 or (name != "migration_rate" and value == 0):
            raise ValueError(f"{name} must be positive (migration_rate may be zero)")
    if "post_decline_fraction" in data and data["post_decline_fraction"] > 1:
        raise ValueError("post_decline_fraction must be at most 1")
    if "post_growth_scaling" in data and data["post_growth_scaling"] < 1:
        raise ValueError("post_growth_scaling must be at least 1")
    if "migration_rate" in data and data["migration_rate"] > 1:
        raise ValueError("migration_rate must be at most 1")
    if (
        challenge == 2
        and data["generations_since_migration"] > data["generations_since_split"]
    ):
        raise ValueError("Contact cannot predate the split")
    return {"parameters": data}


def validate_sweeps(path: Path, challenge: int, dataset: str) -> dict:
    chromosome = "7" if dataset == "testing" else "9"
    intervals = []
    for number, line in enumerate(path.read_text().splitlines(), 1):
        if not line.strip() or line.startswith(("#", "track ", "browser ")):
            continue
        fields = line.split()
        if len(fields) != 3:
            raise ValueError(f"BED line {number}: expected chromosome, start, end")
        contig, start, end = fields
        if contig.removeprefix("chr") != chromosome:
            raise ValueError(
                f"{dataset} sweeps require chromosome {chromosome}, found {contig}"
            )
        start, end = int(start), int(end)
        if start < 0 or end <= start:
            raise ValueError(f"BED line {number}: invalid half-open interval")
        intervals.append((start, end))
    if not intervals or (challenge != 8 and len(intervals) != 1):
        raise ValueError("Expected one sweep interval (one or more for Challenge 8)")
    ordered = sorted(intervals)
    if any(right[0] < left[1] for left, right in zip(ordered, ordered[1:])):
        raise ValueError("Sweep intervals overlap or are duplicated")
    return {"chromosome": chromosome, "intervals": ordered}


def relation_meioses(challenge: int) -> dict[str, int]:
    if challenge == 9:
        return dict(IBD_LABEL_MEIOSES)
    mapping = DOG_LABELS if challenge == 11 else HUMAN_LABELS
    return {
        label: relationship[0] + relationship[1]
        for relationship, label in mapping.items()
    }


def validate_relatedness(path: Path, challenge: int, expected_pairs=None) -> dict:
    with path.open(newline="") as handle:
        sample = handle.read(8192)
        handle.seek(0)
        delimiter = "\t" if "\t" in sample.partition("\n")[0] else ","
        reader = csv.DictReader(handle, delimiter=delimiter)
        headers = reader.fieldnames or []
        allowed = {"indID1", "indID2", "relation", "meioses_count"}
        if (
            len(set(headers)) != len(headers)
            or not {"indID1", "indID2"} <= set(headers)
            or not set(headers) <= allowed
            or not {"relation", "meioses_count"} & set(headers)
        ):
            raise ValueError(
                "Expected unique indID1, indID2, and relation and/or meioses_count headers"
            )
        labels = relation_meioses(challenge)
        pairs, individuals = set(), set()
        for number, row in enumerate(reader, 2):
            if None in row or any(
                value is None or not value.strip() for value in row.values()
            ):
                raise ValueError(f"Row {number}: missing or extra fields")
            left, right = row["indID1"].strip(), row["indID2"].strip()
            if left == right:
                raise ValueError(f"Row {number}: self pair")
            pair = tuple(sorted((left, right)))
            if pair in pairs:
                raise ValueError(f"Row {number}: duplicate unordered pair {pair}")
            pairs.add(pair)
            individuals.update(pair)
            label = row.get("relation")
            if label is not None and label not in labels:
                raise ValueError(
                    f"Row {number}: unaccepted Challenge {challenge} relation {label!r}"
                )
            if "meioses_count" in row:
                try:
                    degree = int(row["meioses_count"])
                except ValueError as error:
                    raise ValueError(
                        f"Row {number}: meioses_count must be an integer"
                    ) from error
                if degree <= 0 or (label is not None and degree != labels[label]):
                    raise ValueError(
                        f"Row {number}: meioses_count conflicts with the relation or is invalid"
                    )
    if not pairs:
        raise ValueError("No relatedness pairs")
    if expected_pairs is not None:
        expected = {tuple(sorted(pair)) for pair in expected_pairs}
        if pairs != expected:
            raise ValueError(
                f"Pair coverage mismatch: {len(expected - pairs)} missing, {len(pairs - expected)} extra"
            )
    return {
        "pairs": len(pairs),
        "individuals": len(individuals),
        "input_pair_coverage_verified": expected_pairs is not None,
    }


def all_pairs(samples: list[str]) -> set[tuple[str, str]]:
    if (
        len(samples) != len(set(samples))
        or len(samples) < 2
        or any(not s.strip() for s in samples)
    ):
        raise ValueError("Expected at least two distinct sample IDs")
    return set(combinations(sorted(samples), 2))


def vcf_samples(path: Path) -> list[str]:
    import pysam

    with pysam.VariantFile(str(path)) as vcf:
        samples = list(vcf.header.samples)
    all_pairs(samples)  # Validate identifiers, including a nonempty cohort.
    return samples


def within_family_pairs(families, samples):
    """Require a partition of the actual cohort, then enumerate its pairs."""
    flattened = [sample for family in families for sample in family]
    if (
        len(flattened) != len(set(flattened))
        or set(flattened) != set(samples)
        or any(not family for family in families)
    ):
        raise ValueError("Families must partition the actual VCF samples exactly once")
    return {
        tuple(sorted(pair)) for family in families for pair in combinations(family, 2)
    }


def validate_submission(
    path: Path, challenge: int, dataset: str, expected_pairs=None
) -> dict:
    if dataset not in {"testing", "final"} or challenge not in range(1, 13):
        raise ValueError("Specify challenge 1–12 and dataset testing or final")
    if challenge <= 4:
        result = validate_demography(path, challenge)
    elif challenge <= 8:
        result = validate_sweeps(path, challenge, dataset)
    else:
        result = validate_relatedness(path, challenge, expected_pairs)
    return {
        "challenge": challenge,
        "dataset": dataset,
        "format_valid": True,
        "scope": "file format and supplied pair coverage; not inference accuracy or input provenance",
        **result,
    }
