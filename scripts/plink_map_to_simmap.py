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
        destination.write("#chr\tpos\tmale_cM\tfemale_cM\n")
        for line in source:
            chromosome, _, cm, position = line.split()[:4]
            destination.write(f"{chromosome}\t{position}\t{cm}\t{cm}\n")


if __name__ == "__main__":
    main()
