from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter
from pathlib import Path

import numpy as np

from ghist_relatedness.labels import HUMAN_LABELS, accepted_label


MEIOSES = {label: up + down for (up, down, _), label in HUMAN_LABELS.items()}


def is_accepted_pedigree_source(source: object) -> bool:
    return "accepted pedigree tuple" in str(source)


def pair_key(row: dict[str, object]) -> tuple[str, str]:
    return tuple(sorted((str(row["indID1"]), str(row["indID2"]))))


def topology_evidence(row: dict[str, object]) -> tuple[str | None, int | None]:
    relationship = row["relationship_tuple"]
    if len(relationship) != 3 or relationship[2] is None:
        return None, None
    up, down, ancestors = relationship
    if not math.isfinite(float(up)) or not math.isfinite(float(down)):
        return None, None
    relationship_tuple = (int(up), int(down), int(ancestors))
    return accepted_label(relationship_tuple, "human"), int(up) + int(down)


def choose_relation(
    probabilities: dict[str, float],
    topology_label: str | None,
    topology_meioses: int | None,
    topology_boost: float | None,
    require_flat_degree_agreement: bool = False,
) -> tuple[str, float, str]:
    if topology_boost is None:
        if topology_label is not None:
            return topology_label, probabilities[topology_label], "accepted pedigree tuple"
        allowed = [
            label
            for label in probabilities
            if topology_meioses is None or MEIOSES[label] == topology_meioses
        ]
        if not allowed:
            allowed = list(probabilities)
        relation = max(allowed, key=probabilities.__getitem__)
        total = sum(probabilities[label] for label in allowed)
        return (
            relation,
            probabilities[relation] / max(total, 1e-12),
            "pair model constrained by pedigree meioses",
        )

    flat_relation = max(probabilities, key=probabilities.__getitem__)
    if require_flat_degree_agreement:
        if (
            topology_label is not None
            and topology_meioses == MEIOSES[flat_relation]
        ):
            return (
                topology_label,
                probabilities[topology_label],
                "background-guarded accepted pedigree tuple",
            )
        return (
            flat_relation,
            probabilities[flat_relation],
            "calibrated pair prediction",
        )

    degree_scores = {
        meioses: sum(
            probability
            for label, probability in probabilities.items()
            if MEIOSES[label] == meioses
        )
        for meioses in set(MEIOSES.values())
    }
    if topology_meioses is not None and topology_meioses in degree_scores:
        degree_scores[topology_meioses] *= topology_boost
    selected_meioses = max(degree_scores, key=degree_scores.__getitem__)
    if topology_label is not None and selected_meioses == topology_meioses:
        relation = topology_label
    else:
        allowed = [
            label for label in probabilities if MEIOSES[label] == selected_meioses
        ]
        relation = max(allowed, key=probabilities.__getitem__)
    source = (
        "background-guarded accepted pedigree tuple"
        if relation == topology_label
        else "background-corrected pair consensus"
    )
    return (
        relation,
        degree_scores[selected_meioses] / max(sum(degree_scores.values()), 1e-12),
        source,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pair-predictions", type=Path, nargs="+", required=True)
    parser.add_argument("--bonsai", type=Path, nargs="+", required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--topology-boost", type=float)
    parser.add_argument("--require-flat-degree-agreement", action="store_true")
    args = parser.parse_args()

    prediction_sets = []
    for path in args.pair_predictions:
        rows = json.loads(path.read_text())
        prediction_sets.append({pair_key(row): row for row in rows})
    pair_keys = sorted(prediction_sets[0])
    if any(set(run) != set(pair_keys) for run in prediction_sets):
        raise ValueError("Sensitivity runs contain different pairs")
    bonsai_runs = []
    for path in args.bonsai:
        payload = json.loads(path.read_text())
        run = {pair_key(row): row for row in payload["relationships"]}
        if not set(run).issubset(pair_keys):
            raise ValueError("Pedigree reconstruction contains unknown pairs")
        bonsai_runs.append(run)

    rows = []
    for pair in pair_keys:
        evidence = [
            topology_evidence(
                run.get(
                    pair,
                    {"relationship_tuple": [float("inf"), float("inf"), None]},
                )
            )
            for run in bonsai_runs
        ]
        topology_label, topology_meioses = evidence[0]
        if topology_label is None or any(item != evidence[0] for item in evidence[1:]):
            topology_label, topology_meioses = None, None
        averaged = {
            label: float(
                np.mean([run[pair]["probabilities"][label] for run in prediction_sets])
            )
            for label in MEIOSES
        }
        relation, confidence, source = choose_relation(
            averaged,
            topology_label,
            topology_meioses,
            args.topology_boost,
            args.require_flat_degree_agreement,
        )
        sensitivity = [
            choose_relation(
                run[pair]["probabilities"],
                topology_label,
                topology_meioses,
                args.topology_boost,
                args.require_flat_degree_agreement,
            )[0]
            for run in prediction_sets
        ]
        rows.append(
            {
                "indID1": pair[0],
                "indID2": pair[1],
                "relation": relation,
                "confidence": confidence,
                "source": source,
                "topology_meioses": topology_meioses,
                "sensitivity_labels": sensitivity,
            }
        )

    by_pair = {pair_key(row): row for row in rows}
    samples = sorted({sample for pair in pair_keys for sample in pair})
    parent = {sample: sample for sample in samples}

    def find(sample: str) -> str:
        while parent[sample] != sample:
            parent[sample] = parent[parent[sample]]
            sample = parent[sample]
        return sample

    def union(left: str, right: str) -> None:
        left_root, right_root = find(left), find(right)
        if left_root != right_root:
            parent[right_root] = left_root

    for row in rows:
        if row["relation"] == "siblings" and is_accepted_pedigree_source(row["source"]):
            union(str(row["indID1"]), str(row["indID2"]))
    components: dict[str, list[str]] = {}
    for sample in samples:
        components.setdefault(find(sample), []).append(sample)
    propagations = 0
    conflicts = 0
    for component in components.values():
        if len(component) < 2:
            continue
        for third in set(samples) - set(component):
            edges = [
                by_pair[pair_key({"indID1": sample, "indID2": third})]
                for sample in component
            ]
            accepted = {
                (str(edge["relation"]), edge["topology_meioses"])
                for edge in edges
                if is_accepted_pedigree_source(edge["source"])
            }
            if len(accepted) > 1:
                conflicts += 1
                continue
            if len(accepted) != 1:
                continue
            relation, meioses = next(iter(accepted))
            for edge in edges:
                if (
                    not is_accepted_pedigree_source(edge["source"])
                    and edge["topology_meioses"] == meioses
                    and MEIOSES[str(edge["relation"])] == meioses
                    and edge["relation"] != relation
                ):
                    edge["relation"] = relation
                    edge["source"] = "full-sibling component propagation"
                    propagations += 1

    stable = [len(set(row["sensitivity_labels"])) == 1 for row in rows]
    degree_stable = [
        len({MEIOSES[label] for label in row["sensitivity_labels"]}) == 1
        for row in rows
    ]
    args.out_dir.mkdir(parents=True, exist_ok=True)
    with (args.out_dir / "candidate.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("indID1", "indID2", "relation"))
        writer.writeheader()
        writer.writerows({key: row[key] for key in writer.fieldnames} for row in rows)
    (args.out_dir / "predictions.json").write_text(json.dumps(rows, indent=2) + "\n")
    summary = {
        "pairs": len(rows),
        "individuals": len(samples),
        "label_counts": dict(Counter(row["relation"] for row in rows)),
        "sensitivity_runs": len(prediction_sets),
        "sensitivity_stability": float(np.mean(stable)),
        "degree_sensitivity_stability": float(np.mean(degree_stable)),
        "unstable_pairs": int(len(stable) - sum(stable)),
        "topology_boost": args.topology_boost,
        "topology_runs": len(bonsai_runs),
        "require_flat_degree_agreement": args.require_flat_degree_agreement,
        "full_sibling_component_propagations": propagations,
        "sibling_component_conflicts": conflicts,
        "source_counts": dict(Counter(row["source"] for row in rows)),
        "submission_count": 0,
    }
    (args.out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
