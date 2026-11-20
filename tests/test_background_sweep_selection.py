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


def test_final_sweep_export_and_exon_filter(tmp_path):
    import pytest
    from scripts.select_background_sweep_candidates import read_exons, write_candidate_bed
    exons = tmp_path / 'exons.bed'
    exons.write_text('chr7\t100\t200\nchr9\t300\t400\n')
    assert read_exons(exons, '9') == [(300, 400)]
    with pytest.raises(ValueError, match='No exons'):
        read_exons(exons, '8')
    with pytest.raises(ValueError, match='No exons'):
        nearest_exonic_position(100, [])
    output = tmp_path / 'final.bed'
    write_candidate_bed(output, [{'exonic_center': 350}], '9', 500)
    assert output.read_text() == 'chr9\t0\t850\n'
