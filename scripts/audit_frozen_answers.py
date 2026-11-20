"""Report public testing accuracy after freezing; never tune or modify estimators."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def audit(root: Path, answers: dict) -> list[dict]:
    reports = []
    for freeze in sorted((root / "experiments").glob("*/task*freeze.json")):
        task = int(freeze.name.split("_")[0][4:])
        if str(task) not in answers:
            continue
        record = json.loads(freeze.read_text())
        prediction = record.get("testing_prediction", record.get("primary"))
        truth = answers[str(task)]
        report = {"challenge": task, "historical_freeze": str(freeze.relative_to(root))}
        if task <= 4:
            renames = (
                {
                    "mainland_population_size": "population1_effective_size",
                    "island_population_size": "population2_effective_size",
                    "generations_since_migration": "generations_since_contact",
                }
                if task == 2
                else {}
            )
            report["parameters"] = {
                name: {
                    "prediction": prediction[renames.get(name, name)],
                    "truth": value,
                    "absolute_relative_error": abs(
                        prediction[renames.get(name, name)] / value - 1
                    ),
                }
                for name, value in truth.items()
            }
        else:
            intervals = prediction.get("regions", [prediction])
            covered = [
                any(
                    int(r["chromosome"]) == chrom
                    and r["start"] <= start
                    and end <= r["end"]
                    for r in intervals
                )
                for chrom, start, end in truth["regions"]
            ]
            report.update(
                {
                    "predicted_intervals": len(intervals),
                    "published_sweep_sites": len(covered),
                    "sites_covered": sum(covered),
                    "coverage_fraction": sum(covered) / len(covered),
                    "note": "Coordinate overlap only; this is not an official competition score",
                }
            )
        reports.append(report)
    return sorted(reports, key=lambda item: item["challenge"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument(
        "--answers",
        type=Path,
        default=Path("benchmarks/holdout/testing_answers_2026_11_20.json"),
    )
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    answers = json.loads(args.answers.read_text())
    report = {
        "protocol": "Post-freeze audit only. Public answers must not enter fitting or model selection.",
        "answer_source": str(args.answers),
        "results": audit(args.root, answers["answers"]),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
