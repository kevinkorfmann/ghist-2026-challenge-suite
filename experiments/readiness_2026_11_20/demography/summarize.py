"""Summarize all outcomes without selecting models from simulation truth."""

import json
from pathlib import Path
import numpy as np

root = Path(__file__).resolve().parent


def records(directory):
    result = []
    for path in sorted((root / directory).glob("task*.json")):
        row = json.loads(path.read_text())
        if "configuration" in row:
            result.append((path, row))
    return result


def summarize(directory, planned):
    rows = records(directory)
    groups = {}
    for path, row in rows:
        c = row["configuration"]
        for model, fit in row["fits"].items():
            key = f"task{c['task']}|point{c['point']}|{c['condition']}|{model}"
            groups.setdefault(key, []).append((row, fit, model))
    result = {
        "planned": planned,
        "final_records": len(rows),
        "completed": sum(row["status"] == "completed" for _, row in rows),
        "partial_or_failed": sum(row["status"] != "completed" for _, row in rows),
        "missing": planned - len(rows),
        "groups": {},
        "growth_convention_discrimination": {},
    }
    for _, row in rows:
        c = row["configuration"]
        if c["task"] == 3 and set(row["fits"]) == {"exponential", "sudden"}:
            matched = c["condition"]
            other = "sudden" if matched == "exponential" else "exponential"
            key = f"point{c['point']}|{matched}"
            result["growth_convention_discrimination"].setdefault(key, []).append(
                {
                    "replicate": row["replicate"],
                    "matched_minus_other_log_likelihood": row["fits"][matched][
                        "multinomial_log_likelihood"
                    ]
                    - row["fits"][other]["multinomial_log_likelihood"],
                    "both_selected_restarts_converged": all(
                        fit.get("best_restart_converged") is True
                        for fit in row["fits"].values()
                    ),
                }
            )
    for key, entries in groups.items():
        parameters = list(entries[0][0]["relative_error"][entries[0][2]])
        errors = np.array(
            [
                [row["relative_error"][model][p] for p in parameters]
                for row, _, model in entries
            ]
        )
        result["groups"][key] = {
            "replicates": [row["replicate"] for row, _, _ in entries],
            "ancestry_seeds": [row["seeds"]["ancestry"] for row, _, _ in entries],
            "median_absolute_relative_error": dict(
                zip(parameters, np.median(abs(errors), axis=0).tolist())
            ),
            "maximum_absolute_relative_error": dict(
                zip(parameters, np.max(abs(errors), axis=0).tolist())
            ),
            "signed_relative_errors": [
                row["relative_error"][model] for row, _, model in entries
            ],
            "numerically_supported": [
                fit["optimizer_diagnostics"]["numerically_supported"]
                for _, fit, _ in entries
            ],
            "best_restart_converged": [
                fit.get("best_restart_converged") for _, fit, _ in entries
            ],
            "successful_restarts": [
                fit["successful_restarts"] for _, fit, _ in entries
            ],
            "failed_restart_counts": [
                len(fit["failed_restarts"]) for _, fit, _ in entries
            ],
        }
    return result


summary = {
    name: summarize(name, count)
    for name, count in [
        ("neutral", 22),
        ("growth_optimizer_status", 8),
        ("neutral_10mb", 28),
    ]
}
previous = root.parents[1] / "scientific_2026_11_20/demography"
# Compare previously completed records, retaining old evidence byte-for-byte.
comparison = {}
new = {path.name: row for path, row in records("neutral")}
for path in sorted((previous / "neutral").glob("task*.json")):
    old = json.loads(path.read_text())
    if path.name in new:
        current = new[path.name]
        comparison[path.name] = {
            model: {
                "maximum_absolute_parameter_relative_difference": max(
                    abs(current["fits"][model][p] / old["fits"][model][p] - 1)
                    for p in old["truth"]
                ),
                "log_likelihood_difference": current["fits"][model][
                    "multinomial_log_likelihood"
                ]
                - old["fits"][model]["multinomial_log_likelihood"],
            }
            for model in old["fits"]
        }
summary["comparison_with_previous_completed_runs"] = comparison
summary["growth_optimizer_status_equivalence"] = {}
for path, current in records("growth_optimizer_status"):
    original = new[path.name]
    summary["growth_optimizer_status_equivalence"][path.name] = {
        model: {
            "maximum_absolute_parameter_relative_difference": max(
                abs(current["fits"][model][p] / original["fits"][model][p] - 1)
                for p in original["truth"]
            ),
            "log_likelihood_difference": current["fits"][model][
                "multinomial_log_likelihood"
            ]
            - original["fits"][model]["multinomial_log_likelihood"],
        }
        for model in original["fits"]
    }
summary["formerly_incomplete"] = {}
old_execution = json.loads((previous / "execution.json").read_text())
for record in old_execution["neutral"]["records"]:
    if record["status"] == "completed":
        continue
    c = record["configuration"]
    name = f"task{c['task']}_point{c['point']}_{c['condition']}_rep{record['replicate']}.json"
    summary["formerly_incomplete"][name] = new.get(name, {}).get("status", "missing")
path = root / "bgs/split_q1.json"
if path.exists():
    row = json.loads(path.read_text())
    summary["unscaled_bgs_split"] = {
        "relative_error": row["relative_error"],
        "best_restart_converged": row["fit"]["best_restart_converged"],
        "converged_restarts": row["fit"]["converged_restarts"],
        "finite_restarts": row["fit"]["successful_restarts"],
        "failed_restarts": row["fit"]["failed_restarts"],
        "spectrum": row["spectrum"],
    }
summary["scope"] = [
    "Completed means the numerical workflow returned and checkpointed both fits; it is not an accuracy pass.",
    "Every dataset is an independent linked chromosome; groups are fixed truth point, condition and full replicate.",
    "Counter-v2 10 Mb runs are additional independent seeds, not retuned estimators.",
    "Neutral runs do not test BGS; one BGS segment does not establish broad recovery or coverage.",
    "Missing or failed records remain in denominators; no panel with missing results is treated as fully evaluated.",
    "Growth optimizer status is a separate replay on identical saved spectra; parameter/likelihood equality must be verified before borrowing its diagnostics.",
]
(root / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(
    {
        name: (summary[name]["completed"], summary[name]["missing"])
        for name in ("neutral", "growth_optimizer_status", "neutral_10mb")
    }
)
