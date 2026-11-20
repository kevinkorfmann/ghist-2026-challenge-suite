"""Observation-only single-sweep candidate generation.

Output intervals are provisional until independent biological validation supports
their coverage. The consensus flag is an abstention diagnostic, not confidence.
"""

from __future__ import annotations

import numpy as np

from .recovery import localize_from_scan
from .scan import VariantMatrix, scan_windows


METHOD = "blind_callable_median_three_v1"


def predict_single_sweep(variants: VariantMatrix, accessibility: np.ndarray) -> dict:
    windows = scan_windows(
        variants, accessibility, window_size=100_000, step_size=20_000
    )
    valid = np.flatnonzero(np.isfinite(windows["composite_score"]))
    if not len(valid):
        raise ValueError("No callable finite composite windows")
    peak = valid[np.argmax(windows["composite_score"][valid])]
    localized = localize_from_scan(variants, windows, accessibility=accessibility)
    estimates = {
        "composite": int(windows["center"][peak]),
        "largest_gap_midpoint": localized["largest_gap_midpoint"],
        "flank_ld_100kb": localized["flank_ld_100kb"],
    }
    values = np.asarray(list(estimates.values()), dtype=float)
    center = int(np.median(values))
    span = float(np.ptp(values))
    return {
        "method": METHOD,
        "status": "provisional candidate; biological accuracy not certified",
        "chromosome": variants.contig,
        "contig_length": variants.contig_length,
        "haplotypes": int(variants.haplotypes.shape[1]),
        "accepted_biallelic_sites": int(len(variants.positions)),
        "skipped_records": variants.skipped_records,
        "coordinate_system": "zero-based half-open BED",
        "localization": localized["localization"],
        "estimates": estimates,
        "center": center,
        "interval": [
            max(0, center - 60_000),
            min(variants.contig_length, center + 60_000),
        ],
        "half_width": 60_000,
        "consensus_span": span,
        "consensus_eligible": span <= 50_000,
        "consensus_threshold": 50_000,
        "submission_count": 0,
    }
