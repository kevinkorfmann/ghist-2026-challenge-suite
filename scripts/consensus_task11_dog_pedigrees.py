from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter
from itertools import combinations
from numbers import Real
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
METHOD_VERSION = "task11-legacy-consensus-audited-v2"
CONFIDENCE_DEFINITION = (
    "Mean uncalibrated pair-model probability of the final label across HBD runs; "
    "not a calibrated probability of correctness."
)


def pair_key(row: dict[str, object]) -> tuple[str, str]:
    samples = (row["indID1"], row["indID2"])
    if any(not isinstance(sample, str) or not sample for sample in samples):
        raise ValueError("Pair identifiers must be nonempty strings")
    if samples[0] == samples[1]:
        raise ValueError(f"Self-pair is not allowed: {samples[0]}")
    return tuple(sorted(samples))


def bonsai_label(row: dict[str, object]) -> tuple[str | None, int | None]:
    relationship = row["relationship_tuple"]
    if not isinstance(relationship, (list, tuple)) or len(relationship) != 3:
        raise ValueError("Bonsai relationship_tuple must have three entries")
    up, down, ancestors = relationship
    # Bonsai uses exactly (+inf, +inf, None) for unrelated individuals.
    # Other absent or nonfinite entries are malformed, not an unconstrained fit.
    if ancestors is None:
        if all(
            not isinstance(value, bool)
            and isinstance(value, Real)
            and value == math.inf
            for value in (up, down)
        ):
            return None, None
        raise ValueError(
            "Unsupported Bonsai unrelated sentinel; expected (+inf, +inf, None)"
        )
    if any(
        isinstance(value, bool) or not isinstance(value, Real) for value in relationship
    ):
        raise ValueError("Bonsai relationship distances and ancestors must be numeric")
    if any(
        not math.isfinite(value) or value < 0 or int(value) != value
        for value in relationship
    ):
        raise ValueError(
            "Finite Bonsai relationship entries must be nonnegative integers"
        )
    canonical = (int(up), int(down), int(ancestors))
    return accepted_label(canonical, "dog"), int(up) + int(down)


def validate_probabilities(
    probabilities: dict[str, float], *, complete: bool = True
) -> None:
    if not isinstance(probabilities, dict) or not probabilities:
        raise ValueError("Pair probabilities must be a nonempty mapping")
    if set(probabilities) - set(MEIOSES):
        raise ValueError("Pair probabilities contain unknown labels")
    if complete and set(probabilities) != set(MEIOSES):
        raise ValueError("Pair probabilities must contain every accepted dog label")
    if any(
        isinstance(value, bool)
        or not isinstance(value, Real)
        or not math.isfinite(value)
        or not 0 <= value <= 1
        for value in probabilities.values()
    ):
        raise ValueError(
            "Pair probabilities must be finite numbers between zero and one"
        )
    if complete and not math.isclose(sum(probabilities.values()), 1.0, abs_tol=1e-6):
        raise ValueError("Pair probabilities must sum to one")


def index_pairs(rows: list[dict], description: str) -> dict:
    indexed = {}
    for row in rows:
        pair = pair_key(row)
        if pair in indexed:
            raise ValueError(f"Duplicate unordered pair in {description}: {pair}")
        indexed[pair] = row
    return indexed


def validate_inputs(prediction_rows: list[list[dict]], bonsai_payload: dict) -> tuple:
    if not prediction_rows:
        raise ValueError("At least one HBD prediction run is required")
    prediction_sets = [
        index_pairs(rows, f"HBD run {i}") for i, rows in enumerate(prediction_rows)
    ]
    pair_keys = sorted(prediction_sets[0])
    if not pair_keys:
        raise ValueError("At least one within-family pair is required")
    if any(set(predictions) != set(pair_keys) for predictions in prediction_sets):
        raise ValueError("HBD sensitivity runs contain different pairs")
    for predictions in prediction_sets:
        for row in predictions.values():
            validate_probabilities(row["probabilities"])

    family_by_sample = {}
    expected_pairs = set()
    family_ids = set()
    for family_index, family in enumerate(bonsai_payload["families"]):
        if family["family"] in family_ids:
            raise ValueError("Duplicate family identifier")
        family_ids.add(family["family"])
        if not family["members"]:
            raise ValueError("Families must not be empty")
        for sample in family["members"]:
            if not isinstance(sample, str) or not sample:
                raise ValueError("Family sample identifiers must be nonempty strings")
            if sample in family_by_sample:
                raise ValueError(f"Duplicate sample in family partition: {sample}")
            family_by_sample[sample] = family_index
        expected_pairs.update(combinations(sorted(family["members"]), 2))
    if set(pair_keys) != expected_pairs:
        raise ValueError(
            "Pair predictions do not exactly cover all within-family pairs"
        )
    bonsai = index_pairs(bonsai_payload["relationships"], "Bonsai relationships")
    if set(bonsai) != expected_pairs:
        raise ValueError(
            "Bonsai relationships do not exactly cover all within-family pairs"
        )
    for row in bonsai.values():
        bonsai_label(row)
    return prediction_sets, pair_keys, bonsai, family_by_sample


def choose_from_probabilities(
    probabilities: dict[str, float], meioses: int | None
) -> tuple[str, float]:
    """Resolve exact ties in the fixed MEIOSES label order, never JSON key order."""
    validate_probabilities(probabilities, complete=False)
    allowed = [
        label
        for label in MEIOSES
        if label in probabilities and (meioses is None or MEIOSES[label] == meioses)
    ]
    if not allowed:
        allowed = [label for label in MEIOSES if label in probabilities]
    total = sum(probabilities[label] for label in allowed)
    if total <= 0:
        raise ValueError("No probability support for topology-constrained labels")
    selected = max(allowed, key=probabilities.__getitem__)
    return selected, probabilities[selected] / total


def build_consensus(
    prediction_rows: list[list[dict]], bonsai_payload: dict
) -> tuple[list[dict], dict]:
    prediction_sets, pair_keys, bonsai, family_by_sample = validate_inputs(
        prediction_rows, bonsai_payload
    )
    labels_by_run: list[dict[tuple[str, str], str]] = []
    raw_labels_by_run = [
        {
            pair: choose_from_probabilities(run[pair]["probabilities"], None)[0]
            for pair in pair_keys
        }
        for run in prediction_sets
    ]
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
                np.mean([run[pair]["probabilities"][label] for run in prediction_sets])
            )
            for label in MEIOSES
        }
        if topology_label is not None:
            relation = topology_label
            source = "accepted pedigree tuple"
        else:
            relation, _ = choose_from_probabilities(averaged, topology_meioses)
            source = "HBD-robust pair model constrained by pedigree meioses"
        consensus_rows.append(
            {
                "indID1": pair[0],
                "indID2": pair[1],
                "relation": relation,
                "confidence": averaged[relation],
                "probabilities": averaged,
                "source": source,
                "topology_meioses": topology_meioses,
                "sensitivity_labels": [run[pair] for run in raw_labels_by_run],
                "topology_conditioned_sensitivity_labels": [
                    run[pair] for run in labels_by_run
                ],
            }
        )

    by_pair = {pair_key(row): row for row in consensus_rows}
    degree_two_refinements = 0
    for row in consensus_rows:
        if row["source"] == "accepted pedigree tuple" or row["topology_meioses"] != 2:
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
        # Preserve the historical refinement rule; raw sensitivity is reporting only.
        if (
            len(set(row["topology_conditioned_sensitivity_labels"])) > 1
            and match_fraction < 0.5
        ):
            row["relation"] = "half_sibling"
            row["source"] = "family relation-vector refinement"
            degree_two_refinements += 1

    sibling_propagations = 0
    sibling_component_conflicts = 0
    for family_index, family in enumerate(bonsai_payload["families"]):
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
            if family_by_sample[str(row["indID1"])] != family_index:
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

    # Every mutation above must report support for the final label, never stale
    # confidence for the replaced label or a renormalized topology-only slice.
    for row in consensus_rows:
        row["confidence"] = row["probabilities"][row["relation"]]
    stable = [len(set(row["sensitivity_labels"])) == 1 for row in consensus_rows]
    conditioned_stable = [
        len(set(row["topology_conditioned_sensitivity_labels"])) == 1
        for row in consensus_rows
    ]
    fallback_stable = [
        is_stable
        for row, is_stable in zip(consensus_rows, stable, strict=True)
        if row["source"] != "accepted pedigree tuple"
    ]
    summary = {
        "method_version": METHOD_VERSION,
        "status": "development only; complete-task validation required",
        "submission_ready": False,
        "confidence_definition": CONFIDENCE_DEFINITION,
        "probability_tie_policy": "Exact ties use the fixed MEIOSES label order; ties are not evidence of accuracy.",
        "hbd_stability_definition": "Pair-model argmax before topology or propagation; stability is not accuracy.",
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
        "fallback_hbd_sensitivity_stability": float(np.mean(fallback_stable))
        if fallback_stable
        else None,
        "fallback_hbd_sensitivity_pairs": len(fallback_stable),
        "topology_conditioned_hbd_sensitivity_stability": float(
            np.mean(conditioned_stable)
        ),
        "unstable_pairs": int(len(stable) - sum(stable)),
        "degree_two_relation_vector_refinements": degree_two_refinements,
        "full_sibling_relation_vector_propagations": sibling_propagations,
        "sibling_component_conflicts": sibling_component_conflicts,
        "final_labels_with_zero_pair_probability": sum(
            row["confidence"] == 0 for row in consensus_rows
        ),
        "submission_count": 0,
    }
    return consensus_rows, summary


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Legacy Task 11 consensus, with audited reporting; not a release gate."
    )
    parser.add_argument("--pair-predictions", type=Path, nargs="+", required=True)
    parser.add_argument("--bonsai", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    consensus_rows, summary = build_consensus(
        [json.loads(path.read_text()) for path in args.pair_predictions],
        json.loads(args.bonsai.read_text()),
    )
    args.out_dir.mkdir(parents=True, exist_ok=True)
    with (args.out_dir / "task11_candidate.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("indID1", "indID2", "relation"))
        writer.writeheader()
        writer.writerows(
            {key: row[key] for key in writer.fieldnames} for row in consensus_rows
        )
    (args.out_dir / "predictions.json").write_text(
        json.dumps(consensus_rows, indent=2, allow_nan=False) + "\n"
    )
    (args.out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
