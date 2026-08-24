from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.cluster import SpectralClustering
from sklearn.metrics import adjusted_rand_score


def sample_number(sample: str) -> int:
    return int(sample.rsplit("_", 1)[1])


def affinity_matrix(
    pairs: pd.DataFrame, samples: list[str], column: str
) -> np.ndarray:
    indices = {sample: index for index, sample in enumerate(samples)}
    affinity = np.zeros((len(samples), len(samples)), dtype=float)
    for row in pairs.itertuples(index=False):
        left = indices[row.indID1]
        right = indices[row.indID2]
        affinity[left, right] = affinity[right, left] = float(getattr(row, column))
    np.fill_diagonal(affinity, float(np.max(affinity)))
    return affinity


def canonical_labels(labels: np.ndarray, samples: list[str]) -> np.ndarray:
    groups = []
    for label in np.unique(labels):
        members = np.flatnonzero(labels == label)
        groups.append((min(sample_number(samples[index]) for index in members), members))
    output = np.empty_like(labels)
    for new_label, (_, members) in enumerate(sorted(groups)):
        output[members] = new_label
    return output


def cluster(affinity: np.ndarray, families: int, seed: int) -> np.ndarray:
    return SpectralClustering(
        n_clusters=families,
        affinity="precomputed",
        assign_labels="kmeans",
        random_state=seed,
    ).fit_predict(affinity)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pairs", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--families", type=int, default=8)
    args = parser.parse_args()

    pairs = pd.read_csv(args.pairs)
    samples = sorted(
        set(pairs["indID1"]) | set(pairs["indID2"]), key=sample_number
    )
    solutions: dict[str, np.ndarray] = {}
    for threshold in (10, 20):
        affinity = affinity_matrix(pairs, samples, f"relatedness_ge_{threshold}")
        for seed in range(8):
            name = f"ge_{threshold}_seed_{seed}"
            solutions[name] = canonical_labels(
                cluster(affinity, args.families, seed), samples
            )
    reference = solutions["ge_10_seed_0"]
    agreement = {
        name: float(adjusted_rand_score(reference, labels))
        for name, labels in solutions.items()
    }
    if min(agreement.values()) < 1.0:
        raise RuntimeError("Family partition is not identical across thresholds and seeds")

    families = []
    for label in range(args.families):
        members = [
            sample for sample, assigned in zip(samples, reference, strict=True) if assigned == label
        ]
        families.append(members)
    summary = {
        "state": "COMPLETED",
        "families": families,
        "family_sizes": [len(members) for members in families],
        "minimum_adjusted_rand": min(agreement.values()),
        "stability_runs": len(agreement),
        "submission_count": 0,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
