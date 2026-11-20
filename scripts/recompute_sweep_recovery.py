"""Recompute legacy simulated localizers without access to the truth position."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import tskit

from ghist_single_sweep.recovery import localize_from_scan
from ghist_single_sweep.scan import VariantMatrix


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", type=Path, required=True)
    parser.add_argument("--trees", type=Path, required=True)
    parser.add_argument("--scan", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    old = json.loads(args.record.read_text())
    ts = tskit.load(args.trees)
    matrix = ts.genotype_matrix()
    positions = np.asarray(
        [site.position - old["left"] for site in ts.sites()], dtype=np.int64
    )
    keep = (
        np.all((matrix == 0) | (matrix == 1), axis=1)
        & (matrix.sum(axis=1) > 0)
        & (matrix.sum(axis=1) < matrix.shape[1])
    )
    variants = VariantMatrix(
        "simulated",
        old["sequence_length"],
        positions[keep],
        matrix[keep].astype(np.uint8),
        int((~keep).sum()),
    )
    table = np.genfromtxt(args.scan, delimiter=",", names=True)
    windows = {name: table[name] for name in table.dtype.names}
    localized = localize_from_scan(variants, windows)
    result = dict(old)
    result["estimates"] = dict(old["estimates"])
    for method in ("flank_ld_100kb", "largest_gap_midpoint"):
        result["estimates"][method] = localized[method] + old["left"]
    result["localization"] = localized["localization"]
    result["absolute_error"] = {
        name: abs(value - old["truth"]["position"])
        for name, value in result["estimates"].items()
    }
    result["previous_oracle_estimates"] = old["estimates"]
    result["reanalysis_of_legacy_simulation"] = True
    result["source_record"] = str(args.record)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
