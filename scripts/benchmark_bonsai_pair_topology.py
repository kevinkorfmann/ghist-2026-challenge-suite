from __future__ import annotations

import argparse
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

from ghist_relatedness.labels import accepted_label


PAIR_PATTERN = re.compile(r"^(.+?)(\d+)_")


def pedigree_name(sample: str) -> tuple[str, int]:
    match = PAIR_PATTERN.match(sample)
    if match is None:
        raise ValueError(f"Cannot parse Ped-sim identifier: {sample}")
    return match.group(1), int(match.group(2))


def read_pairs(path: Path) -> dict[tuple[str, int], list[list[object]]]:
    pairs: dict[tuple[str, int], list[list[object]]] = defaultdict(list)
    with path.open() as handle:
        for line in handle:
            fields = line.split()
            sample1, sample2 = fields[:2]
            label1, replicate1 = pedigree_name(sample1)
            label2, replicate2 = pedigree_name(sample2)
            if (label1, replicate1) != (label2, replicate2):
                raise ValueError("A segment joins different simulated pedigrees")
            pairs[(label1, replicate1)].append(
                [
                    1,
                    2,
                    fields[2],
                    int(fields[3]),
                    int(fields[4]),
                    fields[5] == "IBD2",
                    float(fields[8]),
                ]
            )
    return pairs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--segments", type=Path, required=True)
    parser.add_argument("--threshold-cm", type=float, default=5.0)
    parser.add_argument("--pairs-per-label", type=int, default=20)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    from bonsaitree import bonsai

    bonsai.validate_bio = lambda bio, idx: None
    grouped = read_pairs(args.segments)
    selected = []
    for label in sorted({label for label, _ in grouped}):
        options = sorted(key for key in grouped if key[0] == label)
        selected.extend(options[: args.pairs_per_label])

    rows = []
    for label, replicate in selected:
        segments = [
            segment
            for segment in grouped[(label, replicate)]
            if float(segment[-1]) >= args.threshold_cm
        ]
        if not segments:
            rows.append(
                {
                    "truth": label,
                    "replicate": replicate,
                    "relationship_tuple": None,
                    "prediction": None,
                }
            )
            continue
        total_half = sum(float(segment[-1]) for segment in segments if not segment[-2])
        total_full = sum(float(segment[-1]) for segment in segments if segment[-2])
        result = bonsai.build_pedigree(
            segments,
            [
                {"genotype_id": 1, "age": None, "sex": None},
                {"genotype_id": 2, "age": None, "sex": None},
            ],
            ibd_stats={
                frozenset((1, 2)): {
                    "total_half": total_half,
                    "total_full": total_full,
                    "num_half": sum(not segment[-2] for segment in segments),
                    "max_seg_cm": max((float(segment[-1]) for segment in segments), default=0.0),
                }
            },
            max_add_degree=4,
            max_rel_append_types=3,
            ibd_threshold=0,
            remove_distant_threshold=float("inf"),
        )
        relationship = tuple(result["ped_obj"].rels[1][2])
        predicted = accepted_label(relationship, "human")
        rows.append(
            {
                "truth": label,
                "replicate": replicate,
                "relationship_tuple": relationship,
                "prediction": predicted,
            }
        )

    evaluable = [row for row in rows if row["prediction"] is not None]
    exact = sum(row["truth"] == row["prediction"] for row in evaluable)
    summary = {
        "state": "COMPLETED",
        "threshold_cm": args.threshold_cm,
        "pairs": len(rows),
        "accepted_pairs": len(evaluable),
        "accepted_fraction": len(evaluable) / max(len(rows), 1),
        "exact_accuracy_on_accepted": exact / max(len(evaluable), 1),
        "prediction_counts": dict(Counter(row["prediction"] for row in rows)),
        "rows": rows,
        "submission_count": 0,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({key: value for key, value in summary.items() if key != "rows"}, indent=2))


if __name__ == "__main__":
    main()
