"""Audit every planned paired fit, including failures and censored starts."""
from __future__ import annotations

from collections import Counter
import hashlib
import json
import math
import os
from pathlib import Path
import statistics
import tempfile


HERE = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    manifest = json.loads((HERE / "manifest.json").read_text())
    manifest_hash = digest(HERE / "manifest.json")
    counters = Counter()
    comparisons, missing = [], []
    complete_cases = 0
    for case in manifest["cases"]:
        if digest(HERE / case["input"]) != case["input_sha256"]:
            raise ValueError("Manifest input mismatch")
        directory = HERE / "results" / case["name"]
        completion_path = directory / "completion.json"
        if completion_path.is_file():
            completion = json.loads(completion_path.read_text())
            if completion["manifest_sha256"] != manifest_hash or completion["input_sha256"] != case["input_sha256"]:
                raise ValueError("Completion identity mismatch")
            for name, expected in completion["output_sha256"].items():
                if digest(directory / name) != expected:
                    raise ValueError(f"Completed output mismatch: {case['name']}/{name}")
            complete_cases += 1
        for method in case["methods"]:
            arms = {}
            for arm in manifest["arms"]:
                records, states = [], Counter()
                for restart in range(manifest["restarts"]):
                    path = directory / f"{method}_{arm}_{restart:02d}.json"
                    if not path.is_file():
                        state = "NOT_RECORDED"
                    else:
                        record = json.loads(path.read_text())
                        for field, expected in {
                            "case_name": case["name"], "method": method, "arm": arm, "restart": restart,
                            "manifest_sha256": manifest_hash, "input_sha256": case["input_sha256"],
                            "start_log_parameters": case["log_uniform_starts"][restart],
                        }.items():
                            if record[field] != expected:
                                raise ValueError(f"Fit identity mismatch: {path}, {field}")
                        state = record["state"]
                        if state not in {"RUNNING", "COMPLETED", "FAILED"}:
                            raise ValueError("Unknown fit state")
                        if state == "COMPLETED":
                            if not math.isfinite(record["objective"]):
                                raise ValueError("Completed objective is not finite")
                            if record["optimizer"]["converged"] != (record["optimizer"]["status"] == 0):
                                raise ValueError("Optimizer status disagrees with convergence flag")
                            neighborhood = record["neighborhood"]
                            changes = [probe["objective_change"] for probe in neighborhood["probes"]]
                            if len(changes) != 24 or not all(math.isfinite(value) for value in changes):
                                raise ValueError("Incomplete/nonfinite neighborhood probes")
                            improvement = max(0.0, -min(changes))
                            if improvement != neighborhood["maximum_tested_improvement"]:
                                raise ValueError("Neighborhood improvement does not match raw probes")
                            passes = improvement <= manifest["finite_neighborhood"]["maximum_objective_improvement"]
                            if passes != neighborhood["passes_finite_neighborhood_check"]:
                                raise ValueError("Neighborhood flag does not match raw probes")
                        records.append(record)
                    states[state] += 1
                    counters[state] += 1
                    if state in {"NOT_RECORDED", "RUNNING"}:
                        missing.append({"case": case["name"], "method": method, "arm": arm, "restart": restart, "state": state})
                finite = [row for row in records if row["state"] == "COMPLETED"]
                best = min(finite, key=lambda row: row["objective"]) if finite else None
                arms[arm] = {
                    "planned_restarts": manifest["restarts"], "states": dict(states),
                    "finite_restarts": len(finite),
                    "reported_converged_restarts": sum(row["optimizer"]["converged"] for row in finite),
                    "exceptions": [{"restart": row["restart"], "failure": row["failure"]} for row in records if row["state"] == "FAILED"],
                    "abnormal_or_limited_restarts": [{"restart": row["restart"], "optimizer": row["optimizer"]} for row in finite if not row["optimizer"]["converged"]],
                    "best": None if best is None else {
                        "restart": best["restart"], "parameters": best["parameters"], "objective": best["objective"],
                        "optimizer": best["optimizer"],
                        "restarts_within_0_1": sum(row["objective"] - best["objective"] <= 0.1 for row in finite),
                        "maximum_tested_neighbor_improvement": best["neighborhood"]["maximum_tested_improvement"],
                        "finite_neighborhood_check": best["neighborhood"]["passes_finite_neighborhood_check"],
                        "axis_slopes": best["neighborhood"]["axis_slopes"],
                    },
                }
            fine, coarse = arms["fine"]["best"], arms["coarse"]["best"]
            complete = completion_path.is_file() and all(sum(arms[arm]["states"].get(state, 0) for state in ("COMPLETED", "FAILED")) == manifest["restarts"] for arm in arms)
            criteria = {
                "complete_planned_comparison": complete,
                "no_finite_fit_loss": arms["fine"]["finite_restarts"] >= arms["coarse"]["finite_restarts"],
                "no_material_best_objective_deterioration": fine is not None and coarse is not None and fine["objective"] - coarse["objective"] <= manifest["comparison_criteria"]["maximum_best_objective_deterioration"],
                "selected_fine_best_reports_converged": fine is not None and fine["optimizer"]["converged"],
                "selected_fine_has_restart_agreement": fine is not None and fine["restarts_within_0_1"] >= manifest["comparison_criteria"]["minimum_restarts_within_0_1_of_best"],
                "selected_fine_passes_finite_neighborhood": fine is not None and fine["finite_neighborhood_check"],
            }
            comparisons.append({
                "case": case["name"], "method": method, "arms": arms, "criteria": criteria,
                "all_predeclared_criteria_met": all(criteria.values()),
                "fine_minus_coarse_best_objective": fine["objective"] - coarse["objective"] if fine is not None and coarse is not None else None,
            })
    aggregates = {}
    for group in ("all", "100mb", "10mb", "moments", "dadi"):
        rows = [row for row in comparisons if group == "all" or row["case"].startswith(group) or row["method"] == group]
        differences = [row["fine_minus_coarse_best_objective"] for row in rows if row["fine_minus_coarse_best_objective"] is not None]
        aggregates[group] = {
            "planned_method_problems": len(rows),
            "criteria_met_method_problems": sum(row["all_predeclared_criteria_met"] for row in rows),
            "arms": {arm: {
                "planned_fits": len(rows) * manifest["restarts"],
                "finite_fits": sum(row["arms"][arm]["finite_restarts"] for row in rows),
                "reported_converged_fits": sum(row["arms"][arm]["reported_converged_restarts"] for row in rows),
                "selected_best_converged": sum(bool(row["arms"][arm]["best"] and row["arms"][arm]["best"]["optimizer"]["converged"]) for row in rows),
                "selected_best_passes_neighborhood": sum(bool(row["arms"][arm]["best"] and row["arms"][arm]["best"]["finite_neighborhood_check"]) for row in rows),
                "selected_best_has_restart_agreement": sum(bool(row["arms"][arm]["best"] and row["arms"][arm]["best"]["restarts_within_0_1"] >= 2) for row in rows),
            } for arm in manifest["arms"]},
            "fine_best_improved": sum(value < 0 for value in differences),
            "fine_minus_coarse_best_objective": None if not differences else {
                "minimum": min(differences), "median": statistics.median(differences), "maximum": max(differences)},
        }
    summary = {
        "manifest_sha256": manifest_hash,
        "planned_cases": manifest["planned_cases"], "completed_cases_with_verified_hashes": complete_cases,
        "planned_method_problems": manifest["planned_method_problems"], "planned_fits": manifest["planned_fits"],
        "fit_states": dict(counters), "censored_or_not_recorded": missing,
        "full_length_cases_available_at_freeze": manifest["full_length_cases_available_at_freeze"],
        "full_length_cases_unavailable_at_freeze": manifest["full_length_cases_unavailable_at_freeze"],
        "criteria_met_method_problems": sum(row["all_predeclared_criteria_met"] for row in comparisons),
        "criteria_failure_counts": dict(Counter(key for row in comparisons for key, passed in row["criteria"].items() if not passed)),
        "aggregates": aggregates,
        "counterexamples": [{"case": row["case"], "method": row["method"],
                             "failed_criteria": [key for key, passed in row["criteria"].items() if not passed]}
                            for row in comparisons if not row["all_predeclared_criteria_met"]],
        "comparisons": comparisons,
        "scope": "Frozen available-spectrum numerical comparison; finite-neighborhood probes are not stationarity/KKT certification; failed starts retained",
        "production_default_changed": False, "original_outputs_replaced": False, "truth_errors_used": False,
    }
    with tempfile.NamedTemporaryFile(mode="w", dir=HERE, delete=False) as handle:
        handle.write(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    os.replace(handle.name, HERE / "summary.json")
    print(json.dumps({key: value for key, value in summary.items() if key not in {"comparisons", "censored_or_not_recorded"}}, indent=2))


if __name__ == "__main__":
    main()
