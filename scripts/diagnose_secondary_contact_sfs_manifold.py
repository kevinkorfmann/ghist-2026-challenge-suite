from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path

import numpy as np

from ghist_secondary_contact.spectrum import generator_matched_secondary_contact_spectrum
from ghist_secondary_contact.vcf_fast import scan_secondary_contact_vcf


def _normalized_log_spectrum(values: np.ndarray) -> np.ndarray:
    transformed = np.log1p(np.asarray(values, dtype=np.float64))
    return transformed / (np.linalg.norm(transformed) + 1e-12)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--chunks", required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--reference", type=float, nargs=6, required=True)
    parser.add_argument("--neighbors", type=int, default=25)
    args = parser.parse_args()
    observed_spectrum, metadata = scan_secondary_contact_vcf(args.vcf)
    observed = _normalized_log_spectrum(observed_spectrum).ravel()
    nearest: list[tuple[float, dict]] = []
    all_distances = []
    for path_string in sorted(glob.glob(args.chunks)):
        path = Path(path_string)
        chunk = np.load(path)
        spectra = chunk["joint_sfs"].reshape(len(chunk["truth"]), -1)
        distances = np.linalg.norm(spectra - observed[None], axis=1)
        all_distances.extend(distances.tolist())
        for row in np.argsort(distances)[: args.neighbors]:
            record = {
                "path": str(path),
                "row": int(row),
                "point_id": int(chunk["point_id"][row]),
                "partition": int(chunk["partition"][row]),
                "design_domain": int(chunk["design_domain"][row]),
                "truth": chunk["truth"][row].tolist(),
                "ancestral_size": float(chunk["ancestral_size"][row]),
                "n_variants": int(chunk["n_variants"][row]),
            }
            nearest.append((float(distances[row]), record))
    nearest.sort(key=lambda item: item[0])
    nearest = nearest[: args.neighbors]
    ancestral, size1, size2, split, contact, migration = args.reference
    expected = generator_matched_secondary_contact_spectrum(
        (
            size1 / ancestral,
            size2 / ancestral,
            split / (2 * ancestral),
            contact / split,
            2 * ancestral * migration,
        ),
        (86, 48),
    )
    expected_counts = (
        4 * ancestral * 8.4e-8 * 100_000_000 * np.clip(np.asarray(expected), 0, None)
    )
    reference_distance = float(
        np.linalg.norm(_normalized_log_spectrum(expected_counts).ravel() - observed)
    )
    all_distances_array = np.asarray(all_distances)
    result = {
        "observed_records": metadata["records"],
        "training_rows": len(all_distances),
        "training_distance_quantiles": {
            str(quantile): float(np.quantile(all_distances_array, quantile))
            for quantile in (0.0, 0.001, 0.01, 0.05, 0.5)
        },
        "reference_parameters": args.reference,
        "reference_distance": reference_distance,
        "reference_distance_percentile": float(
            np.mean(all_distances_array <= reference_distance)
        ),
        "nearest_neighbors": [
            {"distance": distance, **record} for distance, record in nearest
        ],
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
