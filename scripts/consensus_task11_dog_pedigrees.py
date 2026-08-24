from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter
from pathlib import Path

import numpy as np

from ghist_relatedness.labels import accepted_label


MEIOSES = {
    "siblings": 2,
    "half_sibling": 2,
    "avuncular": 3,
    "half_avuncular": 3,
    "avuncular_once_removed": 4,
    "full_1st_cousin": 4,
    "half_1st_cousin": 4,
    "avuncular_twice_removed": 5,
}


def pair_key(row: dict[str, object]) -> tuple[str, str]:
    return tuple(sorted((str(row["indID1"]), str(row["indID2"]))))


def bonsai_label(row: dict[str, object]) -> tuple[str | None, int | None]:
    relationship = row["relationship_tuple"]
    if len(relationship) != 3 or relationship[2] is None:
        return None, None
    up, down, ancestors = relationship
    if not math.isfinite(float(up)) or not math.isfinite(float(down)):
        return None, None
    canonical = (int(up), int(down), int(ancestors))
    return accepted_label(canonical, "dog"), int(up) + int(down)


def choose_from_probabilities(
    probabilities: dict[str, float], meioses: int | None
) -> tuple[str, float]:
    allowed = [
        label for label in probabilities if meioses is None or MEIOSES[label] == meioses
    ]
    if not allowed:
        allowed = list(probabilities)
    total = sum(probabilities[label] for label in allowed)
    selected = max(allowed, key=probabilities.__getitem__)
    return selected, probabilities[selected] / max(total, 1e-12)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pair-predictions", type=Path, nargs="+", required=True)
    parser.add_argument("--bonsai", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()

    prediction_sets = []
    for path in args.pair_predictions:
        rows = json.loads(path.read_text())
        prediction_sets.append({pair_key(row): row for row in rows})
    pair_keys = sorted(prediction_sets[0])
    if any(set(predictions) != set(pair_keys) for predictions in prediction_sets):
        raise ValueError("HBD sensitivity runs contain different pairs")

    bonsai_payload = json.loads(args.bonsai.read_text())
    bonsai = {pair_key(row): row for row in bonsai_payload["relationships"]}
    labels_by_run: list[dict[tuple[str, str], str]] = []
    consensus_rows = []
    for predictions in prediction_sets:
        run_labels = {}
        for pair in pair_keys:
            topology_label, topology_meioses = bonsai_label(bonsai[pair])
            if topology_label is not None:
                run_labels[pair] = topology_label
            else:
                run_labels[pair] = choose_from_probabilities(
                    predictions[pair]["probabilities"], topology_meioses
                )[0]
        labels_by_run.append(run_labels)

    for pair in pair_keys:
        topology_label, topology_meioses = bonsai_label(bonsai[pair])
        averaged = {
            label: float(
                np.mean(
                    [run[pair]["probabilities"][label] for run in prediction_sets]
                )
            )
            for label in MEIOSES
        }
        if topology_label is not None:
            relation = topology_label
            confidence = averaged[relation]
            source = "accepted pedigree tuple"
        else:
            relation, confidence = choose_from_probabilities(averaged, topology_meioses)
            source = "HBD-robust pair model constrained by pedigree meioses"
        consensus_rows.append(
            {
                "indID1": pair[0],
                "indID2": pair[1],
                "relation": relation,
                "confidence": confidence,
                "source": source,
                "topology_meioses": topology_meioses,
                "sensitivity_labels": [run[pair] for run in labels_by_run],
            }
        )

    stable = [len(set(row["sensitivity_labels"])) == 1 for row in consensus_rows]
    args.out_dir.mkdir(parents=True, exist_ok=True)
    with (args.out_dir / "task11_candidate.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("indID1", "indID2", "relation"))
        writer.writeheader()
        writer.writerows(
            {key: row[key] for key in writer.fieldnames} for row in consensus_rows
        )
    (args.out_dir / "predictions.json").write_text(
        json.dumps(consensus_rows, indent=2) + "\n"
    )
    summary = {
        "pairs": len(consensus_rows),
        "families": len(bonsai_payload["families"]),
        "label_counts": dict(Counter(row["relation"] for row in consensus_rows)),
        "accepted_topology_pairs": sum(
            row["source"] == "accepted pedigree tuple" for row in consensus_rows
        ),
        "pair_model_fallbacks": sum(
            row["source"] != "accepted pedigree tuple" for row in consensus_rows
        ),
        "hbd_threshold_runs": len(prediction_sets),
        "hbd_sensitivity_stability": float(np.mean(stable)),
        "unstable_pairs": int(len(stable) - sum(stable)),
        "submission_count": 0,
    }
    (args.out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
