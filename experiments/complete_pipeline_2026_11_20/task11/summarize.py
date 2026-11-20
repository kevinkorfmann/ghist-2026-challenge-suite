"""Rebuild Task 11 smoke audit without fitting or modifying any predictions."""

from collections import Counter
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def key(row):
    return tuple(sorted((row["indID1"], row["indID2"])))


def main():
    audit = []
    cohorts = []
    verified_files = 0
    for directory in sorted((ROOT / "results").iterdir()):
        result = json.loads((directory / "result.json").read_text())
        protocol = json.loads((directory / "protocol.json").read_text())
        truth = json.loads((directory / "truth.json").read_text())
        consensus = json.loads((directory / "consensus/predictions.json").read_text())
        forest = [
            {
                key(row): row
                for row in json.loads(
                    (directory / f"pair_predictions_hbd{threshold}.json").read_text()
                )
            }
            for threshold in (2, 5, 10)
        ]
        topology = {
            key(row): row
            for row in json.loads((directory / "bonsai.json").read_text())[
                "relationships"
            ]
        }
        kind_by_sample = {
            sample: kind
            for family, kind in zip(
                truth["families"], truth["family_mating_types"], strict=True
            )
            for sample in family
        }
        for row in consensus:
            pair = key(row)
            probabilities = {
                label: sum(run[pair]["probabilities"][label] for run in forest) / 3
                for label in forest[0][pair]["probabilities"]
            }
            audit.append(
                {
                    "seed": truth["seed"],
                    "indID1": pair[0],
                    "indID2": pair[1],
                    "parent_mating_type": kind_by_sample[pair[0]],
                    "truth": "siblings",
                    "forest_hbd2": forest[0][pair]["relation"],
                    "forest_hbd5": forest[1][pair]["relation"],
                    "forest_hbd10": forest[2][pair]["relation"],
                    "mean_sibling_probability": probabilities["siblings"],
                    "topology_tuple": str(topology[pair]["relationship_tuple"]),
                    "consensus_label": row["relation"],
                    "consensus_source": row["source"],
                    "consensus_reported_confidence": row["confidence"],
                    "mean_probability_of_final_label": probabilities[row["relation"]],
                    "correct": row["relation"] == "siblings",
                }
            )
        cohorts.append(
            {
                "seed": truth["seed"],
                "families": result["family_partition"],
                "consensus": result["consensus"],
                "protocol_sha256": result["protocol_sha256"],
                "model_sha256": protocol["provenance"]["model"]["sha256"],
            }
        )
        for relative, expected in result["output_hashes"].items():
            path = directory / relative
            if path.is_file():
                if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
                    raise ValueError(f"Output hash mismatch: {path}")
                verified_files += 1
    with (ROOT / "pair_before_after.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(audit[0]))
        writer.writeheader()
        writer.writerows(audit)
    wrong = [row for row in audit if not row["correct"]]
    summary = {
        "scope": "held-out whole-family operational smoke only; not a biological validation/release gate",
        "job_id": "10850571",
        "independent_cohorts": len(cohorts),
        "independent_families": sum(
            cohort["families"]["truth_families"] for cohort in cohorts
        ),
        "distinct_pairs": len(audit),
        "forest_correct_by_hbd_threshold": {
            str(threshold): sum(
                row[f"forest_hbd{threshold}"] == row["truth"] for row in audit
            )
            for threshold in (2, 5, 10)
        },
        "consensus_correct": len(audit) - len(wrong),
        "consensus_wrong_by_source": dict(
            Counter(row["consensus_source"] for row in wrong)
        ),
        "wrong_with_zero_reported_confidence": sum(
            row["consensus_reported_confidence"] == 0 for row in wrong
        ),
        "wrong_propagations_with_stale_unit_confidence": sum(
            row["consensus_source"] == "full-sibling component propagation"
            and row["consensus_reported_confidence"] == 1
            for row in wrong
        ),
        "cohorts": cohorts,
        "local_retained_output_hashes_verified": verified_files,
        "source_archive_sha256": "be32c7095bfd47951ad203ad78b5a436efecbd86659e3db77cb44d661ab27ee7",
        "submission_ready": False,
        "submission_count": 0,
    }
    (ROOT / "audit_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(
        json.dumps(
            {key: value for key, value in summary.items() if key != "cohorts"}, indent=2
        )
    )


if __name__ == "__main__":
    main()
