from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path

import numpy as np

from ghist_bottleneck.sfs_baseline import fit_two_epoch_spectrum
from ghist_bottleneck.dadi_reference import fit_two_epoch_dadi_spectrum
try:
    from scripts.evaluate_bottleneck_chromosomes import _metrics, apply_affine, fit_affine
except ModuleNotFoundError:
    from evaluate_bottleneck_chromosomes import _metrics, apply_affine, fit_affine


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def recover_folded_histogram(encoded: np.ndarray, total: int) -> np.ndarray:
    """Invert the normalized log-count encoding used by summary features."""
    encoded = np.asarray(encoded, dtype=float)
    if total == 0:
        return np.zeros_like(encoded, dtype=np.int64)
    if total < 0 or not np.isfinite(encoded).all() or np.any(encoded < 0):
        raise ValueError("invalid folded histogram encoding")
    if not np.any(encoded > 0):
        raise ValueError("a nonempty histogram cannot have an all-zero encoding")
    lower, upper = 0.0, 1.0
    while np.expm1(upper * encoded).sum() < total:
        upper *= 2
    for _ in range(80):
        middle = (lower + upper) / 2
        if np.expm1(middle * encoded).sum() < total:
            lower = middle
        else:
            upper = middle
    raw = np.expm1(((lower + upper) / 2) * encoded)
    counts = np.rint(raw).astype(np.int64)
    difference = int(total - counts.sum())
    if difference > 0:
        order = np.argsort(-(raw - counts))
        for index in np.resize(order, difference):
            counts[index] += 1
    elif difference < 0:
        for index in np.argsort(raw - counts):
            if counts[index] > 0:
                counts[index] -= 1
                difference += 1
                if difference == 0:
                    break
    if counts.sum() != total or np.any(counts < 0):
        raise ValueError("failed to recover integer folded histogram")
    reconstructed = np.log1p(counts)
    reconstructed /= np.linalg.norm(reconstructed)
    if not np.allclose(reconstructed, encoded, atol=2e-5, rtol=2e-5):
        raise ValueError("recovered histogram does not reproduce its encoding")
    return counts


def histogram_to_folded_sfs(histogram: np.ndarray, n_haplotypes: int) -> np.ndarray:
    edges = np.linspace(0, 0.5, len(histogram) + 1)
    spectrum = np.zeros(n_haplotypes + 1, dtype=float)
    used_bins = []
    for count in range(1, n_haplotypes // 2 + 1):
        index = int(np.searchsorted(edges, count / n_haplotypes, side="right") - 1)
        index = min(index, len(histogram) - 1)
        if index in used_bins:
            raise ValueError("histogram bins do not resolve folded allele counts")
        used_bins.append(index)
        spectrum[count] = histogram[index]
    unused = np.ones(len(histogram), dtype=bool)
    unused[used_bins] = False
    if np.any(histogram[unused] != 0):
        raise ValueError("nonzero mass occurs outside resolvable folded counts")
    return spectrum


def _fit_cell(arguments: tuple[int, np.ndarray, int, str]) -> tuple[int, dict]:
    cell, spectrum, restarts, engine = arguments
    fit = fit_two_epoch_spectrum if engine == "moments" else fit_two_epoch_dadi_spectrum
    result = fit(spectrum, restarts=restarts, seed=910_000 + cell, folded=True)
    return cell, result


def _latent_from_result(result: dict) -> np.ndarray:
    fraction = float(result["post_decline_fraction"])
    return np.array(
        [np.log(float(result["generations"])), np.log(fraction / (1 - fraction))]
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", type=Path, required=True)
    parser.add_argument("--split-record", type=Path, required=True)
    parser.add_argument("--observed-fit", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--restarts", type=int, default=12)
    parser.add_argument("--engine", choices=("moments", "dadi"), default="moments")
    args = parser.parse_args()

    data = np.load(args.dataset)
    cell_ids = np.unique(data["cell_id"])
    first = {int(cell): int(np.flatnonzero(data["cell_id"] == cell)[0]) for cell in cell_ids}
    nominal_cells = np.array(
        [cell for cell in cell_ids if int(data["nuisance"][first[int(cell)], 0]) == 0],
        dtype=np.int64,
    )
    n_haplotypes = int(data["n_haplotypes"])
    spectra = {}
    targets = {}
    for cell in nominal_cells:
        indices = np.flatnonzero(data["cell_id"] == cell)
        histogram = np.zeros(64, dtype=np.int64)
        for index in indices:
            histogram += recover_folded_histogram(
                data["sfs"][index, :64], int(data["n_variants"][index])
            )
        spectra[int(cell)] = histogram_to_folded_sfs(histogram, n_haplotypes)
        targets[int(cell)] = np.asarray(data["targets"][indices[0]], dtype=float)

    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        fitted = dict(
            executor.map(
                _fit_cell,
                [(cell, spectra[cell], args.restarts, args.engine) for cell in sorted(spectra)],
            )
        )
    split_record = json.loads(args.split_record.read_text())
    split_cells = {
        name: np.asarray(
            [cell for cell in values if cell in fitted], dtype=np.int64
        )
        for name, values in split_record["split_cells"].items()
    }
    ordered_cells = np.asarray(sorted(fitted), dtype=np.int64)
    cell_index = {int(cell): index for index, cell in enumerate(ordered_cells)}
    prediction = np.asarray([_latent_from_result(fitted[int(cell)]) for cell in ordered_cells])
    target = np.asarray([targets[int(cell)] for cell in ordered_cells])
    indices = {
        name: np.asarray([cell_index[int(cell)] for cell in cells], dtype=np.int64)
        for name, cells in split_cells.items()
    }
    scale = target[indices["calibration"]].std(axis=0) + 1e-6
    matrix, intercept = fit_affine(
        prediction[indices["calibration"]], target[indices["calibration"]]
    )
    calibrated = apply_affine(prediction, matrix, intercept)
    metrics = {}
    for name, values in (("raw", prediction), ("calibrated", calibrated)):
        metrics[name] = {
            split: _metrics(
                values[indices[split]],
                target[indices[split]],
                scale,
                bootstrap_seed=8300 + offset,
            )
            for offset, split in enumerate(("selection", "evaluation"))
        }
    raw_selection = metrics["raw"]["selection"]
    calibrated_selection = metrics["calibrated"]["selection"]
    raw_score = float(np.mean(raw_selection["normalized_rmse"]))
    calibrated_score = float(np.mean(calibrated_selection["normalized_rmse"]))
    natural_guard = (
        calibrated_selection["generation_median_absolute_relative_error"]
        <= 1.05 * raw_selection["generation_median_absolute_relative_error"]
        and calibrated_selection["fraction_median_absolute_error"]
        <= 1.05 * raw_selection["fraction_median_absolute_error"]
    )
    chosen = "calibrated" if calibrated_score < 0.98 * raw_score and natural_guard else "raw"
    observed = json.loads(args.observed_fit.read_text())
    observed_raw = _latent_from_result(observed)
    observed_chosen = (
        apply_affine(observed_raw[None], matrix, intercept)[0]
        if chosen == "calibrated"
        else observed_raw
    )
    observed_fraction = float(1 / (1 + np.exp(-observed_chosen[1])))
    serializable = {
        "dataset": str(args.dataset),
        "dataset_sha256": sha256(args.dataset),
        "split_record": str(args.split_record),
        "split_record_sha256": sha256(args.split_record),
        "observed_fit_sha256": sha256(args.observed_fit),
        "evaluator_sha256": sha256(Path(__file__)),
        "nominal_cells": len(ordered_cells),
        "split_sizes": {name: len(values) for name, values in split_cells.items()},
        "restarts_per_cell": args.restarts,
        "engine": args.engine,
        "affine_matrix": matrix.tolist(),
        "affine_intercept": intercept.tolist(),
        "selection_rule": "calibrated normalized RMSE improves at least 2% and both natural-scale medians remain within 5% of raw",
        "chosen": chosen,
        "metrics": metrics,
        "observed": {
            "source": str(args.observed_fit),
            "raw_generations": int(round(np.exp(observed_raw[0]))),
            "raw_post_decline_fraction": float(1 / (1 + np.exp(-observed_raw[1]))),
            "chosen_generations": int(round(np.exp(observed_chosen[0]))),
            "chosen_post_decline_fraction": observed_fraction,
        },
        "cell_fits": {str(cell): fitted[int(cell)] for cell in ordered_cells},
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(serializable, indent=2) + "\n")
    np.savez_compressed(
        args.out.with_suffix(".npz"),
        cell_id=ordered_cells,
        target_latent=target,
        raw_latent=prediction,
        calibrated_latent=calibrated,
    )
    print(json.dumps({key: value for key, value in serializable.items() if key != "cell_fits"}, indent=2))


if __name__ == "__main__":
    main()
