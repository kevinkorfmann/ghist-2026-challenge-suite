"""Apply the published GHIST 2026 weighted Jaccard formula to retained rows.

This supplements, and does not change, prospectively frozen acceptance gates.
The formula is documented at https://www.synapse.org/Synapse:syn74803798/wiki/640739.
This is an independent implementation, not the organizer's executable scorer.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

from scripts.score_relatedness_cohorts import _index


SOURCE_URL = "https://www.synapse.org/Synapse:syn74803798/wiki/640739"


def _summarize(counts: dict, labels: list[str]) -> dict:
    total = sum(counts[label]["truth"] for label in labels)
    if not total:
        raise ValueError("At least one truth pair is required")
    per_label = {}
    for label in labels:
        row = counts[label]
        denominator = row["true_positive"] + row["false_negative"] + row["false_positive"]
        per_label[label] = {
            **dict(row),
            "jaccard": row["true_positive"] / denominator if denominator else None,
            "weight": row["truth"] / total,
        }
    return {
        "truth_pairs": total,
        "weighted_jaccard": sum(
            row["weight"] * (row["jaccard"] or 0.0) for row in per_label.values()
        ),
        "per_label": per_label,
    }


def evaluate_published_metric(cohorts: list[dict], labels: list[str]) -> dict:
    if (not labels or len(set(labels)) != len(labels)
            or any(not isinstance(label, str) or not label for label in labels)):
        raise ValueError("Supply a nonempty unique label vocabulary")
    if not cohorts:
        raise ValueError("At least one cohort is required")
    identifiers, summaries = set(), []
    pooled = {label: Counter(truth=0, true_positive=0, false_negative=0, false_positive=0)
              for label in labels}
    for cohort in cohorts:
        identifier = cohort["cohort_id"]
        if not isinstance(identifier, str) or not identifier or identifier in identifiers:
            raise ValueError("Cohort identifiers must be nonempty and unique")
        identifiers.add(identifier)
        truth = _index(cohort["truth"], labels, truth=True)
        predictions = _index(cohort["predictions"], labels, truth=False)
        counts = {label: Counter(truth=0, true_positive=0, false_negative=0, false_positive=0)
                  for label in labels}
        for pair, row in truth.items():
            expected = row["relation"]
            inferred = predictions.get(pair, {}).get("relation")
            counts[expected]["truth"] += 1
            counts[expected]["true_positive" if inferred == expected else "false_negative"] += 1
        for pair, row in predictions.items():
            inferred = row["relation"]
            if inferred is not None and truth.get(pair, {}).get("relation") != inferred:
                counts[inferred]["false_positive"] += 1
        missing = len(truth.keys() - predictions.keys())
        extra = len(predictions.keys() - truth.keys())
        abstained = sum(row["relation"] is None for row in predictions.values())
        summaries.append({
            "cohort_id": identifier, **_summarize(counts, labels),
            "missing_pairs": missing, "extra_pairs": extra, "abstained_pairs": abstained,
            "pair_coverage_passed": missing == extra == abstained == 0,
        })
        for label in labels:
            pooled[label].update(counts[label])
    return {
        "schema": "ghist2026-published-relatedness-metric-v1",
        "source_url": SOURCE_URL,
        "formula": "sum_k N_k * TP_k / (TP_k + FN_k + FP_k) / sum_k N_k",
        "scope": "Independent implementation of the published formula; not organizer-executable parity or scientific release authorization",
        "pooled": _summarize(pooled, labels), "cohorts": summaries,
        "pair_coverage_passed": all(row["pair_coverage_passed"] for row in summaries),
        "release_authorized": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    original = args.input.read_bytes()
    data = json.loads(original)
    report = evaluate_published_metric(data["cohorts"], data["labels"])
    if args.input.read_bytes() != original:
        raise RuntimeError("Scoring input changed during evaluation")
    report.update(input_sha256=hashlib.sha256(original).hexdigest(),
                  scorer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  row_validator_sha256=hashlib.sha256(Path(__file__).with_name("score_relatedness_cohorts.py").read_bytes()).hexdigest())
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("x") as handle:
        handle.write(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps(report["pooled"]))


if __name__ == "__main__":
    main()
