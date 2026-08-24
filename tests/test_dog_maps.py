from __future__ import annotations

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


class DogMapTests(unittest.TestCase):
    def test_converts_average_maps(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for chromosome in range(1, 39):
                (root / f"chr{chromosome}_average_canFam3.1.txt").write_text(
                    f"chr\tpos\trate\tcM\n{chromosome}\t10\t1.0\t0\n"
                    f"{chromosome}\t20\t1.0\t1.5\n"
                )
            output = root / "map.txt"
            subprocess.run(
                [
                    sys.executable,
                    "scripts/dog_maps_to_plink.py",
                    "--map-dir",
                    str(root),
                    "--output",
                    str(output),
                ],
                check=True,
            )
            lines = output.read_text().splitlines()
        self.assertEqual(lines[:3], ["1\t.\t0.000000000\t1", "1\t.\t0.000000000\t10", "1\t.\t1.500000000\t20"])
        self.assertEqual(len(lines), 114)


if __name__ == "__main__":
    unittest.main()
