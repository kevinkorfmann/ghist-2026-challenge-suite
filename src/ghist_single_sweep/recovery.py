"""Truth-blind sweep localization and validation provenance.

Simulation truth may be used to score predictions, never to choose the region
in which the predictor searches. Historical recovery records lack this contract
and must be regenerated before they support a new scientific freeze.
"""

from __future__ import annotations

import numpy as np

from .scan import VariantMatrix, accessible_bases, flank_ld_scan, merge_intervals


def localize_from_scan(
    variants: VariantMatrix, windows: dict, radius=750_000, accessibility=None
):
    """Select a search region from the observed scan, without a truth argument."""
    valid = np.flatnonzero(np.isfinite(windows["composite_score"]))
    if not len(valid):
        raise ValueError("No finite composite scan locations")
    index = valid[np.argmax(windows["composite_score"][valid])]
    center = int(windows["center"][index])
    region = [max(0, center - radius), min(variants.contig_length, center + radius)]
    if accessibility is not None:
        accessibility = merge_intervals(accessibility)
        keep = (
            accessible_bases(variants.positions, variants.positions + 1, accessibility)
            > 0
        )
        variants = VariantMatrix(
            variants.contig,
            variants.contig_length,
            variants.positions[keep],
            variants.haplotypes[keep],
            variants.skipped_records,
        )
    positions = variants.positions[
        (variants.positions >= region[0]) & (variants.positions < region[1])
    ]
    if len(positions) < 2:
        raise ValueError("Selected region has fewer than two polymorphic sites")
    gaps = np.diff(positions.astype(np.int64))
    valid_gaps = np.ones(len(gaps), dtype=bool)
    if accessibility is not None:
        # A missing-data interval is not evidence of reduced diversity.
        valid_gaps = (
            accessible_bases(positions[:-1], positions[1:] + 1, accessibility)
            == gaps + 1
        )
    available = np.flatnonzero(valid_gaps)
    if not len(available):
        raise ValueError("No fully callable variant gaps in selected region")
    gap = int(available[np.argmax(gaps[available])])
    candidates = np.arange(region[0], region[1], 1000, dtype=np.int64)
    if accessibility is not None:
        valid_candidates = (
            accessible_bases(candidates - 100_000, candidates, accessibility) >= 80_000
        ) & (
            accessible_bases(candidates, candidates + 100_000, accessibility) >= 80_000
        )
        candidates = candidates[valid_candidates]
    ld = flank_ld_scan(variants, candidates, flank_width=100_000)
    finite = np.flatnonzero(np.isfinite(ld["flank_ld_ratio"]))
    if not len(finite):
        raise ValueError("No finite linkage estimates in scan-selected region")
    best_ld = finite[np.argmax(ld["flank_ld_ratio"][finite])]
    return {
        "largest_gap_midpoint": float(
            (int(positions[gap]) + int(positions[gap + 1])) / 2
        ),
        "flank_ld_100kb": int(ld["candidate_position"][best_ld]),
        "largest_gap_length": int(gaps[gap]),
        "localization": {
            "region_source": "full_scan_composite_maximum",
            "region": region,
            "radius": radius,
            "truth_used_for_inference": False,
        },
    }


def validate_blind_records(records: list[dict]) -> None:
    if not records:
        raise ValueError("Recovery contains no records")
    seeds = [record.get("seed") for record in records]
    if any(seed is None for seed in seeds) or len(set(seeds)) != len(seeds):
        raise ValueError("Recovery must contain distinct simulation seeds")
    for record in records:
        localization = record.get("localization", {})
        if (
            localization.get("truth_used_for_inference") is not False
            or localization.get("region_source") != "full_scan_composite_maximum"
        ):
            raise ValueError(
                "Recovery must be regenerated with truth-blind localization"
            )
