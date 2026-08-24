from __future__ import annotations

import gzip
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


class HapIbdPreparationTests(unittest.TestCase):
    def test_imputes_missing_alleles_with_major_allele(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "input.vcf.gz"
            output = root / "output.vcf.gz"
            with gzip.open(source, "wt") as handle:
                handle.write("##fileformat=VCFv4.2\n")
                handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\n")
                handle.write("1\t1\t.\tA\tG\t.\tPASS\t.\tGT\t1|1\t.|0\n")
            subprocess.run(
                [
                    sys.executable,
                    "scripts/prepare_hapibd_vcf.py",
                    "--input",
                    str(source),
                    "--output",
                    str(output),
                ],
                check=True,
            )
            with gzip.open(output, "rt") as handle:
                record = [line for line in handle if not line.startswith("#")][0]
        self.assertTrue(record.rstrip().endswith("1|1\t1|0"))


if __name__ == "__main__":
    unittest.main()
