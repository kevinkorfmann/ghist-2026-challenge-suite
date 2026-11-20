"""Truth-blind sweep localization and validation provenance.

Simulation truth may be used to score predictions, never to choose the region
in which the predictor searches. Historical recovery records lack this contract
and must be regenerated before they support a new scientific freeze.
"""

from __future__ import annotations

import numpy as np

from .scan import VariantMatrix, flank_ld_scan


def localize_from_scan(variants: VariantMatrix, windows: dict, radius=750_000):
    """Select a search region from the observed scan, without a truth argument."""
    valid = np.flatnonzero(np.isfinite(windows["composite_score"]))
    if not len(valid):
        raise ValueError("No finite composite scan locations")
    index = valid[np.argmax(windows["composite_score"][valid])]
    center = int(windows["center"][index])
    region = [max(0, center - radius), min(variants.contig_length, center + radius)]
    positions = variants.positions[
        (variants.positions >= region[0]) & (variants.positions < region[1])
    ]
    if len(positions) < 2:
        raise ValueError("Selected region has fewer than two polymorphic sites")
    gaps = np.diff(positions.astype(np.int64))
    gap = int(np.argmax(gaps))
    candidates = np.arange(region[0], region[1], 1000, dtype=np.int64)
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
