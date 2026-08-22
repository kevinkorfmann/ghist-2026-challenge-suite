from __future__ import annotations

import unittest

import numpy as np

from ghist_split_migration_bgs.spectrum import project_joint_spectrum


class SplitMigrationSpectrumTests(unittest.TestCase):
    def test_projection_preserves_site_count(self):
        spectrum = np.zeros((9, 7), dtype=np.int64)
        spectrum[1, 2] = 7
        spectrum[5, 3] = 11
        projected = project_joint_spectrum(spectrum, (6, 4))
        self.assertEqual(projected.shape, (7, 5))
        self.assertAlmostEqual(float(projected.sum()), 18.0)

    def test_projection_rejects_oversized_sample(self):
        with self.assertRaises(ValueError):
            project_joint_spectrum(np.zeros((5, 5)), (6, 4))


if __name__ == "__main__":
    unittest.main()
