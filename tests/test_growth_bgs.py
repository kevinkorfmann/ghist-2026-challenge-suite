from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from ghist_growth_bgs.audit import (
    audit_growth_resources,
    coding_distance,
    load_coding_intervals,
    load_recombination_map,
)


class GrowthBackgroundSelectionTests(unittest.TestCase):
    def test_resource_loaders_and_joint_audit(self):
        vcf_text = """##fileformat=VCFv4.2
##contig=<ID=2L,length=100>
##contig=<ID=2R,length=100>
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral allele">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb
2L\t11\t.\tA\tG\t.\tPASS\tAA=A\tGT\t0|1\t0|0
2L\t50\t.\tC\tT\t.\tPASS\tAA=T\tGT\t0|1\t1|1
2R\t10\t.\tG\tA\t.\tPASS\tAA=G\tGT\t0|0\t0|1
"""
        bed_text = "2L\t10\t20\n2L\t20\t25\n2R\t5\t15\n"
        map_text = (
            "Chromosome\tPosition(bp)\tRate(cM/Mb)\n"
            "2L\t0\t0.0\n2L\t40\t1.0\n"
            "2R\t0\t0.5\n2R\t50\t2.0\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            vcf = root / "tiny.vcf"
            bed = root / "coding.bed"
            recombination_map = root / "map.txt"
            vcf.write_text(vcf_text)
            bed.write_text(bed_text)
            recombination_map.write_text(map_text)
            coding = load_coding_intervals(bed)
            recombination = load_recombination_map(recombination_map)
            report = audit_growth_resources(vcf, bed, recombination_map)

        self.assertEqual(coding["2L"], [(10, 25)])
        self.assertEqual(coding_distance(coding["2L"], 10), 0)
        self.assertEqual(coding_distance(coding["2L"], 30), 5)
        self.assertEqual(recombination["2R"][1].tolist(), [0.5, 2.0])
        self.assertEqual(report["samples"], 2)
        self.assertEqual(report["sequence_length"], 200)
        self.assertEqual(report["biallelic_sites"], 3)
        self.assertEqual(report["segregating_biallelic_sites"], 3)
        self.assertEqual(report["ancestral_orientation"], {"reference": 2, "alternate": 1})
        self.assertFalse(report["nonsynonymous_annotation_present"])
        self.assertEqual(report["distance_strata"]["coding_overlap"]["biallelic_sites"], 2)


if __name__ == "__main__":
    unittest.main()
