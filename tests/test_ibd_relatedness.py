from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from ghist_relatedness.ibd import (
    pair_features,
    read_segment_rows,
    read_sex_map,
    relationship_from_sample,
)


class IbdRelatednessTests(unittest.TestCase):
    def test_reads_ped_sim_and_builds_pair_features(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "segments.seg"
            path.write_text(
                "siblings1_g2-b1-i1 siblings1_g2-b2-i1 1 0 10 IBD1 0 10 10\n"
                "siblings1_g2-b1-i1 siblings1_g2-b2-i1 1 10 20 IBD2 10 20 10\n"
            )
            map_path = Path(directory) / "map.txt"
            map_path.write_text("1 0 0 0\n1 20 5 15\n")
            records = pair_features(read_segment_rows(path), 100.0, read_sex_map(map_path))
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]["n_segments"], 2.0)
        self.assertAlmostEqual(records[0]["relatedness"], 0.15)
        self.assertAlmostEqual(records[0]["total_male_map_cm"], 5.0)
        self.assertAlmostEqual(records[0]["total_female_map_cm"], 15.0)
        self.assertEqual(relationship_from_sample(str(records[0]["indID1"])), "siblings")
        self.assertIn("ibd1_bin_0", records[0])
        self.assertIn("ibd1_ordered_19", records[0])
        self.assertIn("chromosome_total_ordered_21", records[0])


if __name__ == "__main__":
    unittest.main()


def test_pair_features_merge_reversed_pair_order():
    rows = [dict(indID1=a, indID2=b, chr='1', type='IBD1', length='10', start_cM='0')
            for a, b in [('a', 'b'), ('b', 'a')]]
    result = pair_features(rows, 100.0)
    assert len(result) == 1
    assert result[0]['total_ibd1_cm'] == 20.0
    assert result[0]['n_segments'] == 2


def test_pair_features_reject_invalid_segments():
    import pytest
    row = dict(indID1='a', indID2='b', chr='1', type='IBD1', length='10', start_cM='0')
    for changes in ({'length': 'nan'}, {'type': 'HBD'}, {'indID2': 'a'}):
        with pytest.raises(ValueError):
            pair_features([{**row, **changes}], 100.0)
    with pytest.raises(ValueError, match='genome_cm'):
        pair_features([row], 0.0)
