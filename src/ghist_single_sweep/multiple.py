"""Variable-count sweep candidates with independent-null calibration.

These functions do not establish power under background selection. Calibration
reports must come from independent, matched null simulations, and recovery
must be evaluated separately before a candidate can be frozen.
"""

from __future__ import annotations

import math

import numpy as np


def calibrate_maximum(null_maxima, alpha=0.05):
    maxima = np.asarray(null_maxima, dtype=float)
    if maxima.ndim != 1 or not np.isfinite(maxima).all() or np.any(maxima < 0):
        raise ValueError("Null maxima must be a finite nonnegative vector")
    if not 0 < alpha < 1:
        raise ValueError("alpha must be between zero and one")
    # Finite-sample order-statistic threshold. Use strict > at detection time
    # so ties cannot inflate the nominal error rate under exchangeability.
    rank = math.ceil((len(maxima) + 1) * (1 - alpha))
    if rank > len(maxima) or not len(maxima):
        raise ValueError("Too few independent null replicates for the requested alpha")
    return float(np.sort(maxima)[rank - 1])


def detect_regions(
    report,
    *,
    threshold,
    separation,
    half_width,
    exons,
    chromosome_length,
    max_exon_snap,
):
    """Select every calibrated peak; never require a preset number of sweeps."""
    from .peaks import (
        nearest_exonic_position,
        separated_peaks,
    )

    report = np.asarray(report, dtype=float)
    if not math.isfinite(threshold) or threshold < 0:
        raise ValueError("Threshold must be finite and nonnegative")
    if (
        any(
            not isinstance(v, (int, np.integer)) or v <= 0
            for v in (separation, half_width, chromosome_length)
        )
        or not isinstance(max_exon_snap, (int, np.integer))
        or max_exon_snap < 0
    ):
        raise ValueError("Invalid localization settings")
    if not exons or any(
        start < 0 or end <= start or end > chromosome_length for start, end in exons
    ):
        raise ValueError("Exons must be valid intervals within the chromosome")
    peaks = []
    for index in separated_peaks(report, separation):
        row = report[index]
        # RAiSD centers are 1-based VCF coordinates; exon BED is 0-based.
        position = int(row[0]) - 1
        if not 0 <= position < chromosome_length:
            raise ValueError("RAiSD coordinate lies outside the chromosome")
        if row[6] <= threshold:
            continue
        center = nearest_exonic_position(position, exons)
        if abs(center - position) > max_exon_snap:
            continue
        peaks.append(
            {
                "center": center,
                "score": float(row[6]),
                "raisd_position": int(row[0]),
                "start": max(0, center - half_width),
                "end": min(chromosome_length, center + half_width + 1),
            }
        )
    intervals = []
    for peak in sorted(peaks, key=lambda p: p["start"]):
        if intervals and peak["start"] <= intervals[-1][1]:
            intervals[-1][1] = max(intervals[-1][1], peak["end"])
        else:
            intervals.append([peak["start"], peak["end"]])
    return {
        "peaks": peaks,
        "intervals": intervals,
        "peak_count": len(peaks),
        "region_count": len(intervals),
    }
