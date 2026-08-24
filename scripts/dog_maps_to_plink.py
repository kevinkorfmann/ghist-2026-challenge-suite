from __future__ import annotations

import argparse
import csv
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--map-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as destination:
        for chromosome in range(1, 39):
            path = args.map_dir / f"chr{chromosome}_average_canFam3.1.txt"
            with path.open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            first_position = int(rows[0]["pos"])
            if first_position > 1:
                destination.write(f"{chromosome}\t.\t0.000000000\t1\n")
            for row in rows:
                destination.write(
                    f"{chromosome}\t.\t{float(row['cM']):.9f}\t{int(row['pos'])}\n"
                )


if __name__ == "__main__":
    main()
