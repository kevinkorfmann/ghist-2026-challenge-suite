from __future__ import annotations

import gzip
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


class HavanaExonTests(unittest.TestCase):
    def test_extracts_and_merges_havana_exons(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            gtf = root / "small.gtf.gz"
            bed = root / "exons.bed"
            summary = root / "summary.json"
            with gzip.open(gtf, "wt") as handle:
                handle.write("7\thavana\texon\t11\t20\t.\t+\t.\tgene_id \"a\";\n")
                handle.write("7\thavana\texon\t18\t30\t.\t+\t.\tgene_id \"b\";\n")
                handle.write("7\tensembl\texon\t40\t50\t.\t+\t.\tgene_id \"c\";\n")
            subprocess.run(
                [
                    sys.executable,
                    "scripts/extract_havana_exons.py",
                    "--gtf",
                    str(gtf),
                    "--bed",
                    str(bed),
                    "--summary",
                    str(summary),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            report = json.loads(summary.read_text())
            bed_text = bed.read_text()
        self.assertEqual(bed_text, "chr7\t10\t30\n")
        self.assertEqual(report["merged_intervals"], 1)


if __name__ == "__main__":
    unittest.main()
