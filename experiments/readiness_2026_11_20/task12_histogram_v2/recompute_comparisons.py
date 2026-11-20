"""Audit frozen Task12 CSVs and recompute only the exact-count seed comparison."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from scripts.fit_human_family_relationships import LABELS, relationship_label
from scripts.task12_histogram_protocol import compare_runs, sha256, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--confirmation", type=Path, required=True)
    parser.add_argument("--csv-dir", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    primary = json.loads(args.confirmation.read_text())
    protocol = json.loads(args.protocol.read_text())
    if primary["protocol_sha256"] != sha256(args.protocol):
        raise ValueError("Primary confirmation/protocol mismatch")
    files, seed_pairs = {}, {}
    for run in primary["runs"]:
        seed = run["seed"]
        index = protocol["confirmation"]["seeds"].index(seed)
        for name, metrics in run["models"].items():
            path = args.csv_dir / f"seed{index}_regime{run['regime_index']}_{name}.csv"
            digest = sha256(path)
            if digest != metrics["csv_sha256"]:
                raise ValueError(f"Frozen CSV changed: {path}")
            confusion = {label: dict.fromkeys(LABELS, 0) for label in LABELS}
            pairs = set()
            with path.open(newline="") as handle:
                reader = csv.DictReader(handle)
                if reader.fieldnames != ["indID1", "indID2", "relation"]:
                    raise ValueError("CSV schema mismatch")
                for row in reader:
                    pair = (row["indID1"], row["indID2"])
                    label = relationship_label(pair[0])
                    if (
                        pair in pairs
                        or pair[0] >= pair[1]
                        or relationship_label(pair[1]) != label
                    ):
                        raise ValueError("CSV pair identity/label mismatch")
                    pairs.add(pair)
                    confusion[label][row["relation"]] += 1
            if (
                confusion != metrics["confusion"]
                or len(pairs) != metrics["segment_present_pairs"]
            ):
                raise ValueError("Frozen CSV/report confusion or count mismatch")
            if seed in seed_pairs and seed_pairs[seed] != pairs:
                raise ValueError("Paired models/regimes do not share the same pairs")
            seed_pairs[seed] = pairs
            files[path.name] = digest
    if len(files) != 60 or len(seed_pairs) != 4:
        raise ValueError("Incomplete four-seed/five-regime/three-model evidence")
    corrected = {
        name: compare_runs(primary["runs"], name, protocol["confirmation"]["gate"])
        for name in primary["comparisons"]
    }
    report = {
        "status": "POST_RUN_SCORER_CORRECTION_NO_MODEL_SELECTION",
        "primary_confirmation_sha256": sha256(args.confirmation),
        "primary_protocol_sha256": sha256(args.protocol),
        "primary_source_sha256": protocol["source_sha256"][
            "scripts/task12_histogram_protocol.py"
        ],
        "corrected_source_sha256": sha256(
            Path(__file__).resolve().parents[3] / "scripts/task12_histogram_protocol.py"
        ),
        "recomputation_script_sha256": sha256(__file__),
        "correction": "Classify each seed gain from exact integer correct-count differences over equal paired/regime denominators; floating zero remains a tie.",
        "primary_comparisons": primary["comparisons"],
        "corrected_comparisons": corrected,
        "aggregate_metrics_unchanged": primary["aggregate"],
        "primary_gate_decisions_unchanged": all(
            corrected[name]["gate"] == primary["comparisons"][name]["gate"]
            for name in corrected
        ),
        "verified_csv_sha256": files,
        "segment_present_pairs_by_seed": {
            str(seed): len(pairs) for seed, pairs in seed_pairs.items()
        },
        "generated_pairs": 28000,
        "segment_present_pairs": sum(map(len, seed_pairs.values())),
        "retraining": False,
        "new_predictions": False,
        "candidate_promoted": False,
        "full_pipeline_validated": False,
    }
    write_json(args.output, report)
    print(
        json.dumps(
            {
                "csv_count": len(files),
                "primary_gate_decisions_unchanged": report[
                    "primary_gate_decisions_unchanged"
                ],
                "report_sha256": sha256(args.output),
            }
        )
    )


if __name__ == "__main__":
    main()
