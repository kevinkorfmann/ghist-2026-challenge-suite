from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from ghist_bottleneck.arg_features import gaussian_tmrca_hazard_features
from ghist_bottleneck.features import FeatureCollection


def _array(translation: np.lib.npyio.NpzFile, *names: str) -> np.ndarray:
    for name in names:
        if name in translation:
            return np.asarray(translation[name])
    raise KeyError(f"translation is missing all expected arrays: {', '.join(names)}")


def _item_first(values: np.ndarray, n_items: int) -> np.ndarray:
    if values.shape[0] == n_items:
        return values
    if values.ndim >= 2 and values.shape[1] == n_items:
        return np.moveaxis(values, 1, 0)
    raise ValueError(
        f"cannot identify item axis in posterior shape {values.shape}; "
        f"expected {n_items} items"
    )


def _translation_blocks(
    translation: np.lib.npyio.NpzFile,
    index_map: np.ndarray,
    collection: FeatureCollection,
) -> tuple[np.ndarray, np.ndarray]:
    n_blocks = int(index_map[:, 0].max(initial=-1)) + 1
    if "block_starts" in translation and "block_ends" in translation:
        starts = np.asarray(translation["block_starts"], dtype=float)
        ends = np.asarray(translation["block_ends"], dtype=float)
        if starts.shape != (n_blocks,) or ends.shape != (n_blocks,):
            raise ValueError("translation block coordinates do not match index_map")
        return starts, ends
    if n_blocks != len(collection.block_starts):
        raise ValueError(
            "translation block coordinates are required when its block count "
            "differs from the feature collection"
        )
    return collection.block_starts, collection.block_ends


def inject(collection_path: Path, translation_path: Path, output: Path) -> None:
    collection = FeatureCollection.load(collection_path)
    translation = np.load(translation_path)
    index_map = np.asarray(translation["index_map"], dtype=int)
    if index_map.ndim != 2 or index_map.shape[1] < 1 or np.any(index_map[:, 0] < 0):
        raise ValueError("index_map must contain nonnegative translation block indices")
    means = _item_first(
        _array(translation, "means", "log_tmrca_mean", "log_mean"), len(index_map)
    )
    variances = _item_first(
        _array(translation, "variances", "log_tmrca_variance", "log_variance"),
        len(index_map),
    )
    if means.shape != variances.shape:
        raise ValueError("posterior means and variances must have the same shape")
    translation_starts, translation_ends = _translation_blocks(
        translation, index_map, collection
    )

    spectra = collection.sfs.copy()
    for block, (start, end) in enumerate(
        zip(collection.block_starts, collection.block_ends, strict=True)
    ):
        overlapping = np.flatnonzero(
            (translation_starts < end) & (translation_ends > start)
        )
        rows = np.flatnonzero(np.isin(index_map[:, 0], overlapping))
        if len(rows):
            spectra[block, -33:-1] = gaussian_tmrca_hazard_features(
                means[rows], variances[rows]
            )
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
