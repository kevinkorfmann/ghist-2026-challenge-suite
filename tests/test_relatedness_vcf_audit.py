from __future__ import annotations

import gzip
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


class RelatednessVcfAuditTests(unittest.TestCase):
    def test_counts_phasing_and_missingness(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            vcf = root / "small.vcf.gz"
            out = root / "audit.json"
            with gzip.open(vcf, "wt") as handle:
                handle.write("##fileformat=VCFv4.2\n")
                handle.write("##contig=<ID=1,length=100>\n")
                handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\n")
                handle.write("1\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\t.|.\n")
            subprocess.run(
                [sys.executable, "scripts/audit_relatedness_vcf.py", "--vcf", str(vcf), "--out", str(out)],
                check=True,
                capture_output=True,
                text=True,
            )
            summary = json.loads(out.read_text())
        self.assertEqual(summary["samples"], 2)
        self.assertEqual(summary["variants"], 1)
        self.assertEqual(summary["phased_nonmissing"], 1)
        self.assertEqual(summary["missing_genotypes"], 1)


if __name__ == "__main__":
    unittest.main()
