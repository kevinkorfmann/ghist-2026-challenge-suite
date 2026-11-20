from __future__ import annotations

import argparse
import json
from pathlib import Path

from ghist_bottleneck.provenance import require_diploid_scale


EXPECTED_CANDIDATES = {
    "baseline",
    "wider",
    "deeper_context",
    "extended_local",
    "compact",
    "mamba_context",
}


def select_candidate(root: Path, *, require_complete: bool = True) -> dict:
    paths = sorted(root.glob("*/experiment.json"))
    records = [json.loads(path.read_text()) for path in paths]
    for record in records:
        require_diploid_scale(record, context="Candidate experiment")
    by_candidate = {record["candidate"]: record for record in records}
    if len(by_candidate) != len(records):
        raise ValueError("candidate names are not unique")
    observed = set(by_candidate)
    if require_complete and observed != EXPECTED_CANDIDATES:
        missing = sorted(EXPECTED_CANDIDATES - observed)
        extra = sorted(observed - EXPECTED_CANDIDATES)
        raise ValueError(f"incomplete screen; missing={missing}, extra={extra}")
    if not records:
        raise ValueError(f"no experiment records below {root}")
    identities = {
        (record["dataset_sha256"], int(record["dataset_bytes"])) for record in records
    }
    if len(identities) != 1:
        raise ValueError("candidate records do not use one identical dataset")
    protocols = {
        json.dumps(record["protocol"], sort_keys=True, separators=(",", ":"))
        for record in records
    }
    if len(protocols) != 1:
        raise ValueError("candidate records do not use one identical protocol")
    ranked = sorted(records, key=lambda record: (float(record["selection_score"]), record["candidate"]))
    winner = ranked[0]
    baseline = by_candidate.get("baseline")
    baseline_score = float(baseline["selection_score"]) if baseline else None
    winner_score = float(winner["selection_score"])
    relative_improvement = (
        (baseline_score - winner_score) / baseline_score
        if baseline_score is not None and baseline_score != 0
        else None
    )
    return {
        "coalescent_ploidy": 2,
        "status": "provisional_validation_winner",
        "winner": winner["candidate"],
        "winner_score": winner_score,
        "baseline_score": baseline_score,
        "relative_improvement_over_baseline": relative_improvement,
        "dataset_sha256": winner["dataset_sha256"],
        "dataset_bytes": int(winner["dataset_bytes"]),
        "protocol": winner["protocol"],
        "ranking": [
            {
                "candidate": record["candidate"],
                "selection_score": float(record["selection_score"]),
                "elapsed_seconds": float(record["elapsed_seconds"]),
                "git_revision": record["git_revision"],
            }
            for record in ranked
        ],
        "promotion_gate": (
            "Fresh simulation shard, five seeds, ablations, and bootstrap intervals are still required."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = select_candidate(args.root)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
