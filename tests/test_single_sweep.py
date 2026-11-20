from __future__ import annotations

import unittest

import numpy as np

from ghist_single_sweep.scan import (
    accessible_bases,
    flank_ld_scan,
    haplotype_homozygosity,
    mean_pairwise_r2,
    scan_windows,
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

    def test_scan_exposes_hard_and_soft_scores(self) -> None:
        from ghist_single_sweep.scan import VariantMatrix

        positions = np.arange(0, 1_000, 20, dtype=np.int32)
        patterns = np.asarray(
            [[0, 0, 0, 1, 1, 1], [0, 0, 1, 0, 1, 1], [0, 1, 0, 1, 0, 1]],
            dtype=np.uint8,
        )
        haplotypes = np.vstack([patterns[index % 3] for index in range(positions.size)])
        variants = VariantMatrix("7", 1_000, positions, haplotypes, 0)
        result = scan_windows(variants, np.asarray([[0, 1_000]]), 200, 100)
        self.assertEqual(result["composite_score"].shape, result["soft_composite_score"].shape)
        self.assertIn("z_high_h2_h1", result)
        self.assertIn("z_conditioned_h2_h1", result)
        inactive = result["z_high_h12"] <= 1.0
        np.testing.assert_array_equal(result["z_conditioned_h2_h1"][inactive], 0.0)
        self.assertTrue(np.isnan(result["composite_score"][-1]))
        self.assertTrue(np.isnan(result["soft_composite_score"][-1]))


if __name__ == "__main__":
    unittest.main()


def test_accessibility_union_is_not_double_counted():
    np.testing.assert_array_equal(accessible_bases(np.array([0]), np.array([30]),
                                  np.array([[0, 20], [10, 30]])), [30])


def test_fixed_sites_do_not_change_tajimas_d():
    assert tajima_d(np.array([1, 2, 3, 0, 6]), 1000, 6) == tajima_d(np.array([1, 2, 3]), 1000, 6)


def test_scan_ignores_variants_outside_accessible_intervals():
    from ghist_single_sweep.scan import VariantMatrix
    positions = np.array([5, 10, 20, 50])
    haplotypes = np.array([[0, 0, 1, 1], [0, 1, 0, 1], [0, 1, 1, 0], [1, 0, 1, 0]], dtype=np.uint8)
    variants = VariantMatrix('9', 100, positions, haplotypes, 0)
    scan = scan_windows(variants, np.array([[10, 40], [60, 100]]), 100, 100)
    assert scan['site_count'][0] == 2
    assert scan['accessible_bases'][0] == 70
    empty = scan_windows(variants, np.empty((0, 2), dtype=int), 100, 100)
    assert empty['site_count'][0] == 0
    assert empty['accessible_bases'][0] == 0


def test_vcf_rejects_legacy_zero_coordinates_and_skips_fixed_sites(tmp_path):
    import pytest
    from ghist_single_sweep.scan import read_phased_biallelic_vcf
    path = tmp_path / 'variants.vcf'
    header = ('##fileformat=VCFv4.2\n##contig=<ID=9,length=100>\n'
              '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
              '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\n')
    path.write_text(header + '9\t0\t.\tA\tC\t.\tPASS\t.\tGT\t0|1\t1|0\n')
    with pytest.raises(ValueError, match='1-based'):
        read_phased_biallelic_vcf(path)
    path.write_text(header + '9\t1\t.\tA\tC\t.\tPASS\t.\tGT\t0|1\t1|0\n'
                    '9\t2\t.\tA\tC\t.\tPASS\t.\tGT\t1|1\t1|1\n')
    result = read_phased_biallelic_vcf(path)
    assert result.positions.tolist() == [0]
    assert result.skipped_records == 1
