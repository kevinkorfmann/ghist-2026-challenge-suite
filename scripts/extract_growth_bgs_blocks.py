from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from ghist_growth_bgs.blocks import extract_block_spectra


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--coding-bed", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--block-size", type=int, default=5_000_000)
    args = parser.parse_args()
    payload = extract_block_spectra(
        args.vcf,
        args.coding_bed,
        block_size=args.block_size,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.out, **payload)
    print(
        f"blocks={payload['spectra'].shape[0]} "
        f"sites={int(payload['accepted_sites'])} "
        f"exposure={int(payload['exposure_bp'].sum())}"
    )


if __name__ == "__main__":
    main()
