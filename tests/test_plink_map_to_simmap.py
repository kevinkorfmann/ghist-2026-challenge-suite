import subprocess
import sys
from pathlib import Path


def test_conversion(tmp_path: Path) -> None:
    source = tmp_path / "map.txt"
    output = tmp_path / "map.simmap"
    source.write_text("1 . 0.5 100\n")
    subprocess.run(
        [
            sys.executable,
            "scripts/plink_map_to_simmap.py",
            "--input",
            str(source),
            "--output",
            str(output),
        ],
        check=True,
    )
    assert output.read_text() == "#chr\tpos\tmale_cM\tfemale_cM\n1\t100\t0.5\t0.5\n"
