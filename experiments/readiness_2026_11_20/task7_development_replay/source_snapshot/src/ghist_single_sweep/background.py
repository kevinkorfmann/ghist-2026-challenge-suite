"""Versioned, observation-only Task 7 inference with explicit callable exposure.

This is a new provisional method, not a reproduction of the archived unmasked
localizer and not a biological calibration. RAiSD must scan the same masked SNPs.
"""

from __future__ import annotations

import numpy as np

from .peaks import highest_peak_position, separated_peaks
from .scan import (
    VariantMatrix,
    accessible_bases,
    flank_ld_scan,
    merge_intervals,
    scan_windows,
)

METHOD = "blind_callable_median_ten_v1"
SETTINGS = {
    "composite_widths": [25_000, 50_000, 100_000, 250_000, 500_000],
    "composite_step_divisor": 5,
    "region_width": 100_000,
    "region_radius": 750_000,
    "grid_step": 1000,
    "flank_widths": [50_000, 100_000, 250_000],
    "minimum_callable_flank_fraction": 0.8,
    "raisd_snp_window": 50,
    "raisd_mask_policy": "entire inclusive SNP span callable",
    "half_width": 80_000,
    "secondary_peak_separation": 100_000,
    "secondary_peak_minimum_distance": 2_000_000,
    "minimum_signal_ratio": 4.0,
}


def masked_variants(variants, accessibility):
    mask = merge_intervals(accessibility)
    if not len(mask) or mask[-1, 1] > variants.contig_length:
        raise ValueError("Accessibility must lie within the declared contig")
    keep = accessible_bases(variants.positions, variants.positions + 1, mask) == 1
    if keep.sum() < 2:
        raise ValueError("Fewer than two callable SNPs")
    return VariantMatrix(
        variants.contig,
        variants.contig_length,
        variants.positions[keep],
        variants.haplotypes[keep],
        variants.skipped_records,
    ), mask


def callable_raisd_windows(report, accessibility, contig_length):
    """Convert RAiSD -R 1-based inclusive start/end to half-open mask spans."""
    report = np.asarray(report, dtype=float)
    if (
        report.ndim != 2
        or not len(report)
        or report.shape[1] != 7
        or not np.isfinite(report).all()
        or np.any(report[:, 3:] < 0)
        or np.any(report[:, :3] != np.floor(report[:, :3]))
        or np.any(report[:, 1] < 1)
        or np.any(report[:, 2] > contig_length)
        or np.any(report[:, 1] > report[:, 0])
        or np.any(report[:, 0] > report[:, 2])
    ):
        raise ValueError("RAiSD needs finite valid seven-column -R window coordinates")
    starts = report[:, 1].astype(np.int64) - 1
    ends = report[:, 2].astype(np.int64)
    keep = accessible_bases(starts, ends, accessibility) == ends - starts
    if not keep.any():
        raise ValueError("No fully callable RAiSD windows")
    return report[keep]


def predict_masked_ten(variants, accessibility, raisd_report):
    variants, mask = masked_variants(variants, accessibility)
    report = callable_raisd_windows(raisd_report, mask, variants.contig_length)
    estimates = {}
    for width in SETTINGS["composite_widths"]:
        scan = scan_windows(
            variants, mask, width, width // SETTINGS["composite_step_divisor"]
        )
        valid = np.flatnonzero(np.isfinite(scan["composite_score"]))
        if not len(valid):
            raise ValueError(f"No finite callable composite windows at width {width}")
        best = valid[np.argmax(scan["composite_score"][valid])]
        estimates[f"composite_{width}"] = int(scan["center"][best])
    anchor = estimates[f"composite_{SETTINGS['region_width']}"]
    region = [
        max(0, anchor - SETTINGS["region_radius"]),
        min(variants.contig_length, anchor + SETTINGS["region_radius"]),
    ]
    local = variants.positions[
        (variants.positions >= region[0]) & (variants.positions < region[1])
    ]
    gaps = np.diff(local.astype(np.int64))
    valid = np.flatnonzero(
        accessible_bases(local[:-1], local[1:] + 1, mask) == gaps + 1
    )
    if not len(valid):
        raise ValueError("No fully callable SNP gaps in scan-selected region")
    gap = valid[np.argmax(gaps[valid])]
    estimates["largest_gap_midpoint"] = float(
        (int(local[gap]) + int(local[gap + 1])) / 2
    )
    for width in SETTINGS["flank_widths"]:
        candidates = np.arange(*region, SETTINGS["grid_step"], dtype=np.int64)
        exposure = SETTINGS["minimum_callable_flank_fraction"] * width
        candidates = candidates[
            (accessible_bases(candidates - width, candidates, mask) >= exposure)
            & (accessible_bases(candidates, candidates + width, mask) >= exposure)
        ]
        ld = flank_ld_scan(variants, candidates, width)
        valid = np.flatnonzero(np.isfinite(ld["flank_ld_ratio"]))
        if not len(valid):
            raise ValueError(f"No finite callable LD candidates at flank width {width}")
        best = valid[np.argmax(ld["flank_ld_ratio"][valid])]
        estimates[f"flank_ld_{width}"] = int(ld["candidate_position"][best])
    estimates["raisd"] = highest_peak_position(report, variants.contig_length)
    values = np.asarray(list(estimates.values()))
    center = int(np.median(values))
    half_width = SETTINGS["half_width"]
    interval = [
        max(0, center - half_width),
        min(variants.contig_length, center + half_width),
    ]
    distant = [
        float(report[i, 6])
        for i in separated_peaks(report, SETTINGS["secondary_peak_separation"])
        if abs(int(report[i, 0]) - 1 - center)
        > SETTINGS["secondary_peak_minimum_distance"]
    ]
    baseline = max(distant) if distant else None
    ratio = (
        float(report[:, 6].max() / baseline)
        if baseline is not None and baseline > 0
        else None
    )
    if ratio is not None and not np.isfinite(ratio):
        ratio = None
    eligible = (
        bool(np.all((values >= interval[0]) & (values < interval[1])))
        and bool(np.ptp(values) <= 2 * half_width)
        and ratio is not None
        and ratio >= SETTINGS["minimum_signal_ratio"]
    )
    return {
        "method": METHOD,
        "settings": SETTINGS,
        "status": "provisional candidate; matched biological recovery absent",
        "chromosome": variants.contig.removeprefix("chr"),
        "contig_length": variants.contig_length,
        "haplotypes": int(variants.haplotypes.shape[1]),
        "callable_sites": int(len(variants.positions)),
        "callable_bases": int(np.sum(mask[:, 1] - mask[:, 0])),
        "localization": {
            "region_source": "full_scan_composite_maximum",
            "region": region,
            "truth_used_for_inference": False,
        },
        "estimates": estimates,
        "center": center,
        "interval": interval,
        "consensus_span": float(np.ptp(values)),
        "consensus_eligible": eligible,
        "raisd_windows_total": len(raisd_report),
        "raisd_windows_callable": len(report),
        "distant_peak_score": baseline,
        "signal_ratio": ratio,
        "coordinate_system": "zero-based half-open BED",
        "submission_count": 0,
    }
