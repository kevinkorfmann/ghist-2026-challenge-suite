from __future__ import annotations

import unittest

import numpy as np

from ghist_single_sweep.scan import (
    accessible_bases,
    flank_ld_scan,
    haplotype_homozygosity,
    mean_pairwise_r2,
    tajima_d,
)


class SingleSweepTests(unittest.TestCase):
    def test_accessible_bases_respects_intersections(self) -> None:
        starts = np.array([0, 10, 20])
        ends = np.array([10, 20, 30])
        intervals = np.array([[5, 25]])
        np.testing.assert_array_equal(accessible_bases(starts, ends, intervals), [5, 10, 5])

    def test_haplotype_homozygosity_detects_dominant_pair(self) -> None:
        haplotypes = np.array(
            [
                [0, 0, 0, 0, 1, 1],
                [0, 0, 0, 0, 1, 1],
                [0, 0, 0, 0, 1, 0],
            ],
            dtype=np.uint8,
        )
        h1, h12, h2_h1 = haplotype_homozygosity(haplotypes)
        self.assertGreater(h12, h1)
        self.assertGreater(h2_h1, 0.0)

    def test_ld_and_tajima_statistics_are_finite(self) -> None:
        haplotypes = np.array(
            [[0, 0, 0, 1, 1, 1], [0, 0, 1, 0, 1, 1], [0, 1, 0, 1, 0, 1]],
            dtype=np.uint8,
        )
        self.assertTrue(np.isfinite(mean_pairwise_r2(haplotypes)))
        self.assertTrue(np.isfinite(tajima_d(haplotypes.sum(axis=1), 1000, 6)))

    def test_flank_ld_scan_returns_candidate_grid(self) -> None:
        from ghist_single_sweep.scan import VariantMatrix

        positions = np.arange(0, 200, 10, dtype=np.int32)
        patterns = np.array(
            [[0, 0, 0, 1, 1, 1], [0, 0, 1, 0, 1, 1]], dtype=np.uint8
        )
        haplotypes = np.vstack([patterns[index % 2] for index in range(positions.size)])
        variants = VariantMatrix("7", 200, positions, haplotypes, 0)
        result = flank_ld_scan(variants, np.array([100]), 100)
        self.assertEqual(result["candidate_position"].tolist(), [100])
        self.assertTrue(np.isfinite(result["flank_ld_ratio"][0]))


if __name__ == "__main__":
    unittest.main()
