from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from ghist_bottleneck.arg_features import tmrca_hazard_features
from ghist_bottleneck.features import FeatureCollection


def item_first(values: np.ndarray, n_items: int) -> np.ndarray:
    if values.shape[0] == n_items:
        return values
    if values.ndim >= 2 and values.shape[1] == n_items:
        return np.moveaxis(values, 1, 0)
    raise ValueError(
        f"cannot identify item axis in TMRCA shape {values.shape}; expected {n_items} items"
    )


def inject(collection_path: Path, translation_path: Path, output: Path) -> None:
    collection = FeatureCollection.load(collection_path)
    translation = np.load(translation_path)
    index_map = np.asarray(translation["index_map"], dtype=int)
    values = item_first(np.asarray(translation["log_tmrca"]), len(index_map))
    spectra = collection.sfs.copy()
    for block in range(len(collection.n_variants)):
        rows = np.flatnonzero(index_map[:, 0] == block)
        if len(rows):
            spectra[block, -33:-1] = tmrca_hazard_features(values[rows])
            spectra[block, -1] = 1.0
    enriched = FeatureCollection(
        windows=collection.windows,
        sfs=spectra,
        n_haplotypes=collection.n_haplotypes,
        sequence_length=collection.sequence_length,
        block_length=collection.block_length,
        block_starts=collection.block_starts,
        block_ends=collection.block_ends,
        n_variants=collection.n_variants,
    )
    enriched.save(output)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--collection", required=True, type=Path)
    parser.add_argument("--translation", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    inject(args.collection, args.translation, args.out)


if __name__ == "__main__":
    main()
