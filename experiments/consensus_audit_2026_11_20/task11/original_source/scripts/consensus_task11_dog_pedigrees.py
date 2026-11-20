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

    by_pair = {pair_key(row): row for row in consensus_rows}
    family_by_sample = {
        sample: family_index
        for family_index, family in enumerate(bonsai_payload["families"])
        for sample in family["members"]
    }
    degree_two_refinements = 0
    for row in consensus_rows:
        if (
            row["source"] == "accepted pedigree tuple"
            or row["topology_meioses"] != 2
        ):
            continue
        family = bonsai_payload["families"][family_by_sample[str(row["indID1"])]]
        matches = []
        for third in family["members"]:
            if third in (row["indID1"], row["indID2"]):
                continue
            left = by_pair[pair_key({"indID1": row["indID1"], "indID2": third})]
            right = by_pair[pair_key({"indID1": row["indID2"], "indID2": third})]
            matches.append(left["relation"] == right["relation"])
        match_fraction = float(np.mean(matches)) if matches else 0.0
        row["sibling_relation_vector_match"] = match_fraction
        if len(set(row["sensitivity_labels"])) > 1 and match_fraction < 0.5:
            row["relation"] = "half_sibling"
            row["source"] = "family relation-vector refinement"
            degree_two_refinements += 1

    sibling_propagations = 0
    sibling_component_conflicts = 0
    for family in bonsai_payload["families"]:
        members = family["members"]
        parent = {sample: sample for sample in members}

        def find(sample: str) -> str:
            while parent[sample] != sample:
                parent[sample] = parent[parent[sample]]
                sample = parent[sample]
            return sample

        def union(left: str, right: str) -> None:
            left_root, right_root = find(left), find(right)
            if left_root != right_root:
                parent[right_root] = left_root

        for row in consensus_rows:
            if family_by_sample[str(row["indID1"])] != family["family"]:
                continue
            if row["relation"] == "siblings" and (
                row["source"] == "accepted pedigree tuple"
                or row.get("sibling_relation_vector_match", 0.0) >= 0.5
            ):
                union(str(row["indID1"]), str(row["indID2"]))
        components: dict[str, list[str]] = {}
        for sample in members:
            components.setdefault(find(sample), []).append(sample)
        for component in components.values():
            if len(component) < 2:
                continue
            for third in set(members) - set(component):
                edges = [
                    by_pair[pair_key({"indID1": sample, "indID2": third})]
                    for sample in component
                ]
                accepted = {
                    (str(edge["relation"]), edge["topology_meioses"])
                    for edge in edges
                    if edge["source"] == "accepted pedigree tuple"
                }
                if len(accepted) > 1:
                    sibling_component_conflicts += 1
                    continue
                if not accepted:
                    continue
                relation, meioses = next(iter(accepted))
                for edge in edges:
                    if (
                        edge["source"] != "accepted pedigree tuple"
                        and edge["topology_meioses"] == meioses
                        and edge["relation"] != relation
                    ):
                        edge["relation"] = relation
                        edge["source"] = "full-sibling component propagation"
                        sibling_propagations += 1

    stable = [len(set(row["sensitivity_labels"])) == 1 for row in consensus_rows]
    fallback_stable = [
        is_stable
        for row, is_stable in zip(consensus_rows, stable, strict=True)
        if row["source"] != "accepted pedigree tuple"
    ]
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
        "fallback_hbd_sensitivity_stability": float(np.mean(fallback_stable)),
        "unstable_pairs": int(len(stable) - sum(stable)),
        "degree_two_relation_vector_refinements": degree_two_refinements,
        "full_sibling_relation_vector_propagations": sibling_propagations,
        "sibling_component_conflicts": sibling_component_conflicts,
        "submission_count": 0,
    }
    (args.out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
