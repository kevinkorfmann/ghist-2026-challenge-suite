"""Opt-in development policy; never exports a submission or changes defaults.

Pair probabilities make a proposal. Raw HBD disagreement, an exact probability
tie, or missing/conflicting topology causes abstention. Topology never replaces
the proposal or propagates labels. This conservative guard is unvalidated.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

if __package__:
    from scripts.consensus_task11_dog_pedigrees import (
        CONFIDENCE_DEFINITION,
        MEIOSES,
        bonsai_label,
        validate_inputs,
    )
else:
    from consensus_task11_dog_pedigrees import (
        CONFIDENCE_DEFINITION,
        MEIOSES,
        bonsai_label,
        validate_inputs,
    )


METHOD_VERSION = "task11-guarded-pair-proposal-development-v1"
SOURCE_FILES = (
    "scripts/consensus_task11_guarded_candidate.py",
    "scripts/consensus_task11_dog_pedigrees.py",
    "src/ghist_relatedness/labels.py",
)


def source_hashes() -> dict[str, str]:
    root = Path(__file__).resolve().parents[1]
    return {
        relative: hashlib.sha256((root / relative).read_bytes()).hexdigest()
        for relative in SOURCE_FILES
    }


def unique_argmax(probabilities: dict[str, float]) -> str | None:
    highest = max(probabilities.values())
    winners = [label for label, value in probabilities.items() if value == highest]
    return winners[0] if len(winners) == 1 else None


def build_guarded_candidate(
    prediction_rows: list[list[dict]], bonsai_payload: dict
) -> tuple[list[dict], dict]:
    runs, pairs, topology, _ = validate_inputs(prediction_rows, bonsai_payload)
    rows = []
    for pair in pairs:
        averaged = {
            label: sum(run[pair]["probabilities"][label] for run in runs) / len(runs)
            for label in MEIOSES
        }
        proposal = unique_argmax(averaged)
        raw_labels = [unique_argmax(run[pair]["probabilities"]) for run in runs]
        topology_label, topology_meioses = bonsai_label(topology[pair])
        reasons = []
        if proposal is None or None in raw_labels:
            reasons.append("pair_probability_tie")
        if len(runs) < 2:
            reasons.append("hbd_sensitivity_not_assessed")
        if len(set(raw_labels)) != 1:
            reasons.append("raw_hbd_disagreement")
        if topology_label is None:
            reasons.append("unmapped_topology")
        elif topology_label != proposal:
            reasons.append("topology_disagreement")
        rows.append(
            {
                "indID1": pair[0],
                "indID2": pair[1],
                "proposed_relation": proposal,
                "proposed_relation_pair_probability": averaged[proposal]
                if proposal
                else None,
                "probabilities": averaged,
                "raw_sensitivity_labels": raw_labels,
                "topology_label": topology_label,
                "topology_meioses": topology_meioses,
                "relation": proposal if not reasons else None,
                "decision": "provisional" if not reasons else "abstain",
                "abstention_reasons": reasons,
            }
        )
    provisional = sum(row["decision"] == "provisional" for row in rows)
    summary = {
        "method_version": METHOD_VERSION,
        "status": "DEVELOPMENT; no independent complete-task confirmation",
        "submission_ready": False,
        "submission_count": 0,
        "submission_csv_exported": False,
        "pairs": len(rows),
        "families": len(bonsai_payload["families"]),
        "hbd_threshold_runs": len(runs),
        "provisional_pairs": provisional,
        "abstained_pairs": len(rows) - provisional,
        "provisional_coverage": provisional / len(rows),
        "abstention_reasons": dict(
            Counter(reason for row in rows for reason in row["abstention_reasons"])
        ),
        "confidence_definition": CONFIDENCE_DEFINITION,
        "policy": "Unique pair-model argmax must agree across HBD runs and with a mapped topology label; otherwise abstain. No topology overrides, constraints or propagation.",
        "limitations": "Proposals and provisional decisions are uncalibrated and unvalidated. Agreement is not proof of correctness; abstention cannot fill required submission pairs.",
    }
    return rows, summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--development-policy", choices=[METHOD_VERSION], required=True)
    parser.add_argument("--pair-predictions", type=Path, nargs="+", required=True)
    parser.add_argument("--bonsai", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    rows, summary = build_guarded_candidate(
        [json.loads(path.read_text()) for path in args.pair_predictions],
        json.loads(args.bonsai.read_text()),
    )
    summary["input_sha256"] = {
        str(path): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in [*args.pair_predictions, args.bonsai]
    }
    summary["source_sha256"] = source_hashes()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "development_predictions.json").write_text(
        json.dumps(rows, indent=2, allow_nan=False) + "\n"
    )
    (args.out_dir / "development_summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n"
    )
    print(json.dumps(summary, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
