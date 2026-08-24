from __future__ import annotations

import gzip
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


class HapIbdSummaryTests(unittest.TestCase):
    def test_summarizes_pairs_and_hbd(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            ibd = root / "segments.ibd.gz"
            hbd = root / "segments.hbd.gz"
            map_path = root / "map.txt"
            out = root / "out"
            with gzip.open(ibd, "wt") as handle:
                handle.write("a\t1\tb\t1\t1\t1\t10\t20\n")
                handle.write("a\t2\tb\t2\t1\t20\t30\t10\n")
            with gzip.open(hbd, "wt") as handle:
                handle.write("a\t1\ta\t2\t1\t1\t10\t5\n")
            map_path.write_text("1\t.\t0\t1\n1\t.\t100\t100\n")
            subprocess.run(
                [
                    sys.executable,
                    "scripts/summarize_hapibd.py",
                    "--ibd",
                    str(ibd),
                    "--hbd",
                    str(hbd),
                    "--map",
                    str(map_path),
                    "--out-dir",
                    str(out),
                ],
                check=True,
                capture_output=True,
            )
            summary = json.loads((out / "summary.json").read_text())
        self.assertEqual(summary["pairs_with_ibd"], 1)
        self.assertEqual(summary["hbd_individuals"], 1)
        self.assertEqual(summary["component_sizes"], [2])


if __name__ == "__main__":
    unittest.main()
