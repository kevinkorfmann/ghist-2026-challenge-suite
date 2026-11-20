"""Versioned Task 12 feature correction; historical models retain v1 semantics.

V2 corrects full length-bin counts rather than treating the 20 longest stored
segments as the complete segment list. The existing marginal count/length
quantile corrections remain estimates, not a reconstructed segment list.
No training or prediction entry point opts into v2 automatically.
"""

from __future__ import annotations

import copy
from collections.abc import Mapping
import math

from scripts.fit_human_background_corrected import (
    correct_records as legacy_correct_records,
    corrected_feature_names,
)
from scripts.fit_human_family_relationships import (
    BIN_EDGES,
    THRESHOLDS,
    background_correct,
)


LEGACY = "task12-top20-v1"
FULL_HISTOGRAM = "task12-full-histogram-v2"
METADATA_KEY = "task12_background_correction"


def finite_nonnegative(value, name: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be finite and nonnegative") from exc
    if isinstance(value, bool) or not math.isfinite(result) or result < 0:
        raise ValueError(f"{name} must be finite and nonnegative")
    return result


def correction_metadata(version: str, *, minimum_cm: int, genome_cm: float) -> dict:
    if version not in (LEGACY, FULL_HISTOGRAM):
        raise ValueError(f"Unsupported Task12 correction version: {version!r}")
    if isinstance(minimum_cm, bool) or minimum_cm not in THRESHOLDS:
        raise ValueError("minimum_cm must be an available segment threshold")
    genome = finite_nonnegative(genome_cm, "genome_cm")
    if genome == 0:
        raise ValueError("genome_cm must be positive")
    if version == FULL_HISTOGRAM and minimum_cm != 5:
        raise ValueError("Task12 full-histogram v2 supports only minimum_cm=5")
    return {"version": version, "minimum_cm": minimum_cm, "genome_cm": genome}


def subtract_shortest_histogram(
    counts, count_offset: float, *, minimum_cm: int = 5
) -> list[float]:
    """Remove count mass from complete eligible bins, shortest bins first.

    Fractional offsets are retained, matching the quantile count correction.
    Bins below the threshold are untouched; excess removal saturates at zero.
    This does not estimate the exact lengths of removed segments within a bin.
    """
    if minimum_cm not in THRESHOLDS or isinstance(minimum_cm, bool):
        raise ValueError("minimum_cm must align with an available histogram threshold")
    if len(counts) != len(BIN_EDGES) - 1:
        raise ValueError("Full histogram must contain every length bin")
    output = [finite_nonnegative(value, "histogram count") for value in counts]
    remaining = finite_nonnegative(count_offset, "count_offset")
    for index, lower in enumerate(BIN_EDGES[:-1]):
        if lower >= minimum_cm:
            removed = min(output[index], remaining)
            output[index] -= removed
            remaining -= removed
    return output


def _validate_raw_records(records, genome_cm: float, role: str) -> None:
    if not records:
        raise ValueError(f"Task12 {role} records must be nonempty")
    for row in records:
        try:
            counts = [
                finite_nonnegative(row[f"length_bin_{i}"], f"{role} histogram count")
                for i in range(len(BIN_EDGES) - 1)
            ]
            for threshold in THRESHOLDS:
                count = finite_nonnegative(
                    row[f"segments_ge_{threshold}"], f"{role} segment count"
                )
                total = finite_nonnegative(
                    row[f"total_cm_ge_{threshold}"], f"{role} total cM"
                )
                relatedness = finite_nonnegative(
                    row[f"relatedness_ge_{threshold}"], f"{role} relatedness"
                )
                expected_count = sum(
                    n
                    for n, lower in zip(counts, BIN_EDGES[:-1], strict=True)
                    if lower >= threshold
                )
                if not math.isclose(count, expected_count, rel_tol=1e-9, abs_tol=1e-8):
                    raise ValueError(
                        f"Task12 {role} raw histogram and segment counts disagree at {threshold} cM"
                    )
                if not math.isclose(
                    total / (2 * genome_cm), relatedness, rel_tol=1e-9, abs_tol=1e-10
                ):
                    raise ValueError(
                        f"Task12 {role} total cM and relatedness disagree at {threshold} cM"
                    )
                if total + 1e-8 < threshold * count or (count == 0 and total != 0):
                    raise ValueError(
                        f"Task12 {role} segment lengths violate the {threshold} cM threshold"
                    )
        except KeyError as exc:
            raise ValueError(
                f"Task12 {role} requires complete raw histogram/count/length features"
            ) from exc


def correct_records_v2(records, reference, *, minimum_cm: int, genome_cm: float):
    """Create v2 features for a future explicitly compatible Task12 model."""
    metadata = correction_metadata(
        FULL_HISTOGRAM, minimum_cm=minimum_cm, genome_cm=genome_cm
    )
    genome = metadata["genome_cm"]
    _validate_raw_records(records, genome, "observed")
    _validate_raw_records(reference, genome, "reference")
    corrected = copy.deepcopy(records)
    offsets = background_correct(corrected, reference, 0.25, genome)
    count_offset = offsets[f"segments_ge_{minimum_cm}"]
    for row in corrected:
        counts = subtract_shortest_histogram(
            [row[f"length_bin_{i}"] for i in range(len(BIN_EDGES) - 1)],
            count_offset,
            minimum_cm=minimum_cm,
        )
        row.update({f"length_bin_{i}": value for i, value in enumerate(counts)})
    # max_cm and ordered lengths remain raw and MUST NOT be model features in v2.
    return corrected


def correct_for_model(
    records, reference, model: Mapping, *, minimum_cm: int, genome_cm: float
):
    """Route by artifact metadata; absent metadata means historical v1 only.

    A v2 artifact must be fitted with exactly the v2 feature contract. A legacy
    baseline blend is deliberately unsupported; it has different semantics.
    """
    if METADATA_KEY not in model:
        return legacy_correct_records(records, reference, minimum_cm, genome_cm)
    metadata = model[METADATA_KEY]
    if not isinstance(metadata, Mapping) or set(metadata) != {
        "version",
        "minimum_cm",
        "genome_cm",
    }:
        raise ValueError("Malformed Task12 correction metadata")
    expected = correction_metadata(
        metadata["version"], minimum_cm=minimum_cm, genome_cm=genome_cm
    )
    if dict(metadata) != expected:
        raise ValueError(
            "Task12 model correction metadata does not match the requested input contract"
        )
    if metadata["version"] == LEGACY:
        return legacy_correct_records(records, reference, minimum_cm, genome_cm)
    if model.get("feature_names") != corrected_feature_names(5):
        raise ValueError(
            "Task12 v2 model must use exactly the 5 cM histogram/count/length feature contract"
        )
    if model.get("corrected_model_weight") != 1.0 or isinstance(
        model.get("corrected_model_weight"), bool
    ):
        raise ValueError(
            "Task12 v2 does not support blending an incompatible historical baseline"
        )
    return correct_records_v2(
        records, reference, minimum_cm=minimum_cm, genome_cm=genome_cm
    )
