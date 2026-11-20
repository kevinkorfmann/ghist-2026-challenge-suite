"""Rebuild the relatedness component summary from the adjacent raw JSON files."""

import hashlib
import json
from pathlib import Path
import statistics


folder = Path(__file__).resolve().parent
groups = {}
raw_results = {}
for path in sorted(folder.glob("*.json")):
    data = json.loads(path.read_text())
    if not data.get("component_only") or "runs" not in data:
        continue
    raw_results[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
    version = (
        "current_contract_candidate" if "current" in path.stem else "historical_frozen"
    )
    for condition, run in data["runs"].items():
        key = f"{data['kind']}:{version}:{data['scenario']}:{condition}"
        groups.setdefault(key, []).append((path.name, data, run))

metrics = [
    "accuracy",
    "balanced_accuracy",
    "meioses_accuracy",
    "ece_10_equal_width_bins",
    "multiclass_brier_sum",
    "negative_log_likelihood",
]
report = {
    "design": (
        "Fixed models and fixed metrics on independent simulation seeds; no public "
        "answers; repeated map/background settings are not independent replicates"
    ),
    "component_only": True,
    "candidate_promoted": False,
    "groups": {},
}
for key, rows in groups.items():
    seeds = [row[1]["seed"] for row in rows]
    if len(seeds) != len(set(seeds)):
        raise ValueError(f"Repeated simulation seed in group {key}")
    output = {
        "independent_seeds": len(rows),
        "source_files": [row[0] for row in rows],
        "seeds": seeds,
    }
    for metric in metrics:
        values = [
            row[2].get("accuracy_counting_probability_failures_wrong", row[2]["accuracy"])
            if metric == "accuracy"
            else row[2][metric]
            for row in rows
        ]
        output[metric] = {
            "mean": statistics.mean(values),
            "minimum": min(values),
            "maximum": max(values),
        }
    output["zero_accepted_probability_pairs"] = sum(
        row[2].get("zero_accepted_probability_pairs", 0) for row in rows
    )
    output["genotyped_pair_pipeline_validated"] = False
    output["retained_pairs_before_probability_rejection"] = sum(
        row[2].get("eligible_pairs_before_probability_failure", row[2]["pairs"])
        for row in rows
    )
    output["expected_pairs_including_unobserved_zero_ibd"] = sum(
        row[1]["expected_pairs_before_zero_ibd_censoring"] for row in rows
    )
    output["class_recall_mean"] = {
        label: statistics.mean(row[2]["recall_by_label"][label] for row in rows)
        for label in rows[0][2]["recall_by_label"]
    }
    report["groups"][key] = output

report["result_sha256"] = raw_results
(folder / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
