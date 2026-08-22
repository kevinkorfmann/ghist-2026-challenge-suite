from __future__ import annotations

import argparse
import json
from pathlib import Path

from ghist_secondary_contact.vcf_fast import scan_secondary_contact_vcf


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    spectrum, population = scan_secondary_contact_vcf(args.vcf)
    result = {
        "population": population,
        "joint_sfs_shape": list(spectrum.shape),
        "joint_sfs_sum": int(spectrum.sum()),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
