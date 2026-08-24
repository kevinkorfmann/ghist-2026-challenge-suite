from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.input.open() as source, args.output.open("w") as destination:
        for line in source:
            if not line.strip() or line.startswith("#"):
                continue
            chromosome, position, male_cm, female_cm = line.split()[:4]
            average_cm = (float(male_cm) + float(female_cm)) / 2.0
            destination.write(f"{chromosome}\t.\t{average_cm:.9f}\t{position}\n")


if __name__ == "__main__":
    main()
