from __future__ import annotations

import unittest

import numpy as np

from scripts.select_background_sweep_candidates import nearest_exonic_position, separated_peaks


class BackgroundSweepSelectionTests(unittest.TestCase):
    def test_peak_separation_and_exon_snapping(self) -> None:
        report = np.asarray(
            [[100, 0, 0, 0, 0, 0, 5], [150, 0, 0, 0, 0, 0, 8], [1_000, 0, 0, 0, 0, 0, 6]],
            dtype=float,
        )
        selected = separated_peaks(report, 500)
        self.assertEqual([int(report[index, 0]) for index in selected], [150, 1_000])
        self.assertEqual(nearest_exonic_position(250, [(100, 200), (300, 400)]), 300)


if __name__ == "__main__":
    unittest.main()
