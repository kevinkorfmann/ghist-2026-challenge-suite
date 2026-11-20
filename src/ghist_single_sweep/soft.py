"""SNP-delimited H12 development candidates, not a calibrated sweep detector.

Garud et al. (2015), doi:10.1371/journal.pgen.1005004, motivate SNP-delimited
windows to avoid changing SNP counts in fixed physical haplotype windows.
The window and allele-frequency settings here require independent validation.
"""

from __future__ import annotations

import numpy as np

from .scan import VariantMatrix, accessible_bases, haplotype_homozygosity


def snp_h12_peak(
    variants: VariantMatrix,
    accessibility: np.ndarray,
    window_snps: int,
    minimum_minor_frequency: float = 0,
) -> dict:
    if window_snps < 2 or not 0 <= minimum_minor_frequency < 0.5:
        raise ValueError("Invalid SNP-window settings")
    counts = variants.haplotypes.sum(axis=1)
    n = variants.haplotypes.shape[1]
    keep = (
        accessible_bases(variants.positions, variants.positions + 1, accessibility) > 0
    ) & (np.minimum(counts, n - counts) >= minimum_minor_frequency * n)
    positions, matrix = variants.positions[keep], variants.haplotypes[keep]
    step = max(10, window_snps // 4)
    best = None
    evaluated = 0
    ties = 0
    for left in range(0, len(positions) - window_snps + 1, step):
        right = left + window_snps
        start, end = int(positions[left]), int(positions[right - 1]) + 1
        if end - start > 500_000:
            continue
        # Never form a haplotype window across inaccessible sequence.
        if (
            accessible_bases(np.asarray([start]), np.asarray([end]), accessibility)[0]
            != end - start
        ):
            continue
        h1, h12, h2_h1 = haplotype_homozygosity(matrix[left:right])
        evaluated += 1
        row = dict(
            position=(start + end) // 2,
            start=start,
            end=end,
            h1=h1,
            h12=h12,
            h2_h1=h2_h1,
        )
        if best is None or h12 > best["h12"]:
            best, ties = row, 1
        elif h12 == best["h12"]:
            ties += 1
    if best is None:
        raise ValueError("No fully callable SNP-delimited windows")
    return {
        **best,
        "window_snps": window_snps,
        "minimum_minor_frequency": minimum_minor_frequency,
        "evaluated_windows": evaluated,
        "maximum_ties": ties,
        "tie_rule": "first genomic window; all maximum ties counted",
        "informative": best["h12"] > 1 / n + 2 / n**2 + 1e-12,
        "truth_used_for_inference": False,
    }
