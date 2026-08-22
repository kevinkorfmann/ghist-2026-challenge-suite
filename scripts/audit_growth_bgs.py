from __future__ import annotations

import argparse
import json
from pathlib import Path

from ghist_growth_bgs.audit import audit_growth_resources


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--coding-bed", type=Path, required=True)
    parser.add_argument("--recombination-map", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = audit_growth_resources(args.vcf, args.coding_bed, args.recombination_map)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
