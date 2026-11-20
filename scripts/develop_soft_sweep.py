"""Evaluate explicitly repurposed DEVELOPMENT simulations; never certify them."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import tskit

from ghist_single_sweep.scan import VariantMatrix
from ghist_single_sweep.soft import snp_h12_peak


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trees", type=Path, required=True)
    parser.add_argument("--record", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise ValueError("Refusing to overwrite development evidence")
    record = json.loads(args.record.read_text())
    ts = tskit.load(args.trees)
    matrix = ts.genotype_matrix()
    positions = np.asarray(
        [s.position - record["left"] for s in ts.sites()], dtype=np.int64
    )
    counts = matrix.sum(axis=1)
    keep = (
        np.all((matrix == 0) | (matrix == 1), axis=1)
        & (counts > 0)
        & (counts < ts.num_samples)
    )
    variants = VariantMatrix(
        "7",
        record["sequence_length"],
        positions[keep],
        matrix[keep].astype(np.uint8),
        int((~keep).sum()),
    )
    intervals = np.asarray([[0, variants.contig_length]])
    predictions = {}
    for window in (50, 100, 200, 400):
        for frequency in (0, 0.05):
            name = f"h12_snps{window}_maf{frequency}"
            predictions[name] = snp_h12_peak(variants, intervals, window, frequency)
    truth = record["truth"]["position"] - record["left"]
    with args.trees.open("rb") as handle:
        tree_digest = hashlib.file_digest(handle, "sha256").hexdigest()
    output = {
        "usage": "development only; all 48 former diagnostics may influence method choice; not held out",
        "seed": record["seed"],
        "case": args.record.stem,
        "generator_scale": record["slim_scaling_factor"],
        "coordinate_origin": record["left"],
        "sequence_length": record["sequence_length"],
        "source_record_sha256": hashlib.sha256(args.record.read_bytes()).hexdigest(),
        "source_tree_sha256": tree_digest,
        "coordinate_contract": "predictions are segment-relative; truth.position is chromosome-absolute",
        "truth": record["truth"],
        "predictions": predictions,
        "absolute_errors": {
            k: abs(v["position"] - truth) for k, v in predictions.items()
        },
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({"case": args.record.stem, "state": "COMPLETED"}))


if __name__ == "__main__":
    main()
