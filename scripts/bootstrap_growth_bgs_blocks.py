from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from ghist_growth_bgs.spectrum import fit_growth_spectrum


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--blocks", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--replicate-start", type=int, required=True)
    parser.add_argument("--replicates", type=int, default=4)
    parser.add_argument("--projection", type=int, default=84)
    parser.add_argument("--restarts", type=int, default=12)
    parser.add_argument("--maxiter", type=int, default=400)
    parser.add_argument("--mutation-rate", type=float, default=3.36e-8)
    parser.add_argument("--seed", type=int, default=75194000)
    args = parser.parse_args()

    archive = np.load(args.blocks)
    spectra = archive["spectra"]
    exposure = archive["exposure_bp"]
    contigs = archive["contig"]
    starts = archive["start"]
    ends = archive["end"]
    block_size = int(archive["block_size"])
    usable = (ends - starts) == block_size
    groups = {
        contig: np.flatnonzero((contigs == contig) & usable)
        for contig in np.unique(contigs)
    }
    if any(len(indices) < 2 for indices in groups.values()):
        raise ValueError("each chromosome arm needs at least two full blocks")

    results = []
    for replicate in range(args.replicate_start, args.replicate_start + args.replicates):
        rng = np.random.default_rng(args.seed + replicate)
        selected = np.concatenate(
            [rng.choice(indices, size=len(indices), replace=True) for indices in groups.values()]
        )
        result = fit_growth_spectrum(
            spectra[selected].sum(axis=0),
            sequence_length=float(exposure[selected].sum()),
            mutation_rate=args.mutation_rate,
            projection=args.projection,
            restarts=args.restarts,
            workers=1,
            maxiter=args.maxiter,
            seed=args.seed + 1009 * replicate,
        )
        result["replicate"] = replicate
        result["selected_blocks"] = selected.tolist()
        results.append(result)
        print(
            f"replicate={replicate} generations={result['generations']:.6f} "
            f"scaling={result['post_growth_scaling']:.9f}",
            flush=True,
        )
    payload = {
        "source": str(args.blocks),
        "projection": args.projection,
        "block_size": block_size,
        "full_blocks_by_contig": {
            str(contig): len(indices) for contig, indices in groups.items()
        },
        "terminal_partial_blocks_excluded_from_resampling": True,
        "results": results,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n")


if __name__ == "__main__":
    main()
