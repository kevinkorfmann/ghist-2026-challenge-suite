"""Regenerate compact summary from immutable independent diagnostic records."""

import json
from pathlib import Path
import numpy as np

root = Path(__file__).resolve().parent
records = [
    json.loads(path.read_text()) for path in sorted((root / "neutral").glob("*.json"))
]
groups = {}
for record in records:
    config = record["configuration"]
    for model, errors in record["relative_error"].items():
        key = f"task{config['task']}|point{config['point']}|{config['condition']}|{model}"
        groups.setdefault(key, []).append(record)
summary = {
    "protocol": "Independent msprime-generated linked 2 Mb datasets; truth-independent fixed starts for Task 2/4, globally random growth starts; two replicate seeds per configuration; no challenge data/answers",
    "planned_runs": 22,
    "completed_runs": len(records),
    "incomplete_runs": 22 - len(records),
    "execution_record": "execution.json",
    "groups": {},
}
for key, rows in groups.items():
    model = key.split("|")[-1]
    parameters = list(rows[0]["relative_error"][model])
    errors = np.array(
        [[row["relative_error"][model][key] for key in parameters] for row in rows]
    )
    summary["groups"][key] = {
        "replicates": len(rows),
        "replicate_ids": [row["replicate"] for row in rows],
        "simulation_seeds": [row["seed"] for row in rows],
        "parameters": parameters,
        "median_absolute_relative_error": dict(
            zip(parameters, np.median(abs(errors), axis=0).tolist())
        ),
        "maximum_absolute_relative_error": dict(
            zip(parameters, np.max(abs(errors), axis=0).tolist())
        ),
        "all_optimizer_diagnostics_supported": all(
            row["fits"][model]["optimizer_diagnostics"]["numerically_supported"]
            for row in rows
        ),
    }
summary["bgs"] = {}
for path in sorted((root / "bgs").glob("*_q*.json")):
    row = json.loads(path.read_text())
    if "relative_error" not in row:
        continue  # Execution-state records are preserved, not successful fits.
    summary["bgs"][path.stem] = {
        "relative_error": row["relative_error"],
        "generator": row["generator"],
        "spectrum": row["spectrum"],
    }
summary["limitations"] = [
    "Diagnostic panel, not a full-prior validation or accuracy certificate",
    "Bounded incomplete runs are preserved in execution.json; completed records alone are a censored subset and cannot estimate overall accuracy",
    "Neutral panel has no BGS and uses infinite-sites mutations",
    "Replicate unit is the complete linked simulated chromosome; no random window splitting",
    "No calibrated uncertainty intervals or coverage estimate from two replicates",
    "BGS checks use one 1 Mb segment and independent seeds; cross-Q differences confound stochasticity and approximation",
    "Two-N BGS burn-in and small segment/sample size require longer-run and whole-chromosome stress evaluation",
    "No current final input was used and no prediction was promoted",
]
(root / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(f"Summarized {len(records)} neutral runs and {len(summary['bgs'])} BGS checks")
