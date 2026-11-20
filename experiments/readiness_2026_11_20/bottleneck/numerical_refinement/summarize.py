"""Audit refinement while selecting finite objectives and counting distinct ancestors."""
from __future__ import annotations

from collections import Counter
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile


HERE = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit_finite(record, config):
    if record["state"] != "COMPLETED" or not math.isfinite(record["objective"]):
        raise ValueError("Expected completed finite candidate")
    if record["optimizer"]["converged"] != (record["optimizer"]["status"] == 0):
        raise ValueError("Optimizer status/flag mismatch")
    neighborhood = record["neighborhood"]
    changes = [probe["objective_change"] for probe in neighborhood["probes"]]
    if len(changes) != 24 or not all(math.isfinite(value) for value in changes):
        raise ValueError("Incomplete or nonfinite neighborhood")
    improvement = max(0.0, -min(changes))
    if improvement != neighborhood["maximum_tested_improvement"]:
        raise ValueError("Improvement differs from raw neighborhood probes")
    if neighborhood["passes_finite_neighborhood_check"] != (improvement <= config["finite_neighborhood"]["maximum_objective_improvement"]):
        raise ValueError("Neighborhood pass flag mismatch")


def select_candidates(original, refined, restarts):
    """Exact ties retain an original candidate; status never influences selection."""
    candidates = original + refined
    best_original = min(original, key=lambda candidate: candidate["objective"])
    best = min(candidates, key=lambda candidate: candidate["objective"])
    by_origin = {
        restart: min(candidate["objective"] for candidate in candidates if candidate["restart"] == restart)
        for restart in range(restarts)
    }
    agreement = sum(value - best["objective"] <= 0.1 for value in by_origin.values())
    return best_original, best, agreement, by_origin


def describe(record):
    return {"arm": record["arm"], "restart": record["restart"], "parameters": record["parameters"],
            "objective": record["objective"], "optimizer": record["optimizer"],
            "maximum_tested_neighbor_improvement": record["neighborhood"]["maximum_tested_improvement"],
            "finite_neighborhood_check": record["neighborhood"]["passes_finite_neighborhood_check"],
            "axis_slopes": record["neighborhood"]["axis_slopes"]}


def main():
    config = json.loads((HERE / "manifest.json").read_text())
    manifest_hash = digest(HERE / "manifest.json")
    for name, expected in config["source_sha256"].items():
        if digest(HERE / name) != expected:
            raise ValueError("Frozen diagnostic source mismatch")
    comparisons, incomplete, counters = [], [], Counter()
    completed_cases = 0
    for case in config["cases"]:
        if digest(HERE / case["input"]) != case["input_sha256"]:
            raise ValueError("Input hash mismatch")
        originals = HERE / "original_results" / case["name"]
        output = HERE / "results" / case["name"]
        for name, expected in case["original_output_sha256"].items():
            if digest(originals / name) != expected:
                raise ValueError("Preserved original endpoint mismatch")
        completion_path = output / "completion.json"
        if completion_path.is_file():
            completion = json.loads(completion_path.read_text())
            if completion["manifest_sha256"] != manifest_hash or completion["input_sha256"] != case["input_sha256"]:
                raise ValueError("Completion identity mismatch")
            for name, expected in completion["output_sha256"].items():
                if digest(output / name) != expected:
                    raise ValueError("Completed output hash mismatch")
            completed_cases += 1
        for method in case["methods"]:
            original = []
            for restart in range(config["restarts"]):
                for arm in ("coarse", "fine"):
                    record = json.loads((originals / f"{method}_{arm}_{restart:02d}.json").read_text())
                    for field, expected in {"case_name": case["name"], "method": method, "arm": arm, "restart": restart,
                                            "manifest_sha256": config["parent_manifest_sha256"],
                                            "input_sha256": case["input_sha256"],
                                            "start_log_parameters": case["log_uniform_starts"][restart]}.items():
                        if record[field] != expected:
                            raise ValueError("Original endpoint identity mismatch")
                    audit_finite(record, config)
                    original.append(record)
            refined, failures, states = [], [], Counter()
            for restart, start in enumerate(case["refinement_starts"][method]):
                path = output / f"{method}_refined_{restart:02d}.json"
                state = "NOT_RECORDED"
                if path.is_file():
                    record = json.loads(path.read_text())
                    for field, expected in {"case_name": case["name"], "method": method, "arm": "refined", "restart": restart,
                                            "manifest_sha256": manifest_hash, "parent_manifest_sha256": config["parent_manifest_sha256"],
                                            "input_sha256": case["input_sha256"], "start_log_parameters": start["log_parameters"],
                                            "original_random_start_log_parameters": case["log_uniform_starts"][restart],
                                            "coarse_endpoint_sha256": start["coarse_endpoint_sha256"]}.items():
                        if record[field] != expected:
                            raise ValueError(f"Refinement identity mismatch: {path}, {field}")
                    state = record["state"]
                    if state not in {"RUNNING", "COMPLETED", "FAILED"}:
                        raise ValueError("Unknown refinement state")
                    if state == "COMPLETED":
                        audit_finite(record, config)
                        coarse = next(candidate for candidate in original if candidate["arm"] == "coarse" and candidate["restart"] == restart)
                        if abs(record["initial_objective"] - coarse["objective"]) > 1e-8:
                            raise ValueError("Coarse initial objective changed")
                        refined.append(record)
                    elif state == "FAILED":
                        failures.append({"restart": restart, "failure": record["failure"]})
                states[state] += 1
                counters[state] += 1
                if state in {"NOT_RECORDED", "RUNNING"}:
                    incomplete.append({"case": case["name"], "method": method, "restart": restart, "state": state})
            best_original, best, agreement, origins = select_candidates(original, refined, config["restarts"])
            original_agreement = sum(min(row["objective"] for row in original if row["restart"] == restart) - best_original["objective"] <= 0.1 for restart in range(config["restarts"]))
            criteria = {
                "all_planned_refinements_complete": completion_path.is_file() and states["COMPLETED"] == config["restarts"],
                "no_objective_degradation_vs_original_union": best["objective"] <= best_original["objective"],
                "selected_best_reports_converged": best["optimizer"]["converged"],
                "selected_best_has_distinct_origin_agreement": agreement >= config["criteria"]["minimum_distinct_origins_within_0_1"],
                "selected_best_passes_finite_neighborhood": best["neighborhood"]["passes_finite_neighborhood_check"],
            }
            comparisons.append({"case": case["name"], "method": method, "refinement_states": dict(states),
                                "exceptions": failures,
                                "refined_converged": sum(row["optimizer"]["converged"] for row in refined),
                                "refined_abnormal_or_limited": [{"restart": row["restart"], "optimizer": row["optimizer"]} for row in refined if not row["optimizer"]["converged"]],
                                "refinement_worsened_ancestor_objective": sum(row["objective"] > row["initial_objective"] for row in refined),
                                "original_union_best": describe(best_original), "selected_best": describe(best),
                                "selected_minus_original_union_objective": best["objective"] - best_original["objective"],
                                "original_distinct_origins_within_0_1": original_agreement,
                                "distinct_origins_within_0_1": agreement,
                                "best_objective_by_original_start": origins,
                                "criteria": criteria, "all_predeclared_criteria_met": all(criteria.values())})
    aggregates = {}
    for group in ("all", "100mb", "10mb", "moments", "dadi"):
        rows = [row for row in comparisons if group == "all" or row["case"].startswith(group) or row["method"] == group]
        aggregates[group] = {"planned_method_problems": len(rows),
                             "criteria_met": sum(row["all_predeclared_criteria_met"] for row in rows),
                             "selected_converged": sum(row["selected_best"]["optimizer"]["converged"] for row in rows),
                             "selected_passes_neighborhood": sum(row["selected_best"]["finite_neighborhood_check"] for row in rows),
                             "selected_has_distinct_origin_agreement": sum(row["distinct_origins_within_0_1"] >= 2 for row in rows),
                             "selected_arm_counts": dict(Counter(row["selected_best"]["arm"] for row in rows)),
                             "objective_improved_vs_original_union": sum(row["selected_minus_original_union_objective"] < 0 for row in rows),
                             "refined_reported_converged": sum(row["refined_converged"] for row in rows)}
    summary = {"manifest_sha256": manifest_hash, "parent_manifest_sha256": config["parent_manifest_sha256"],
               "planned_cases": config["planned_cases"], "completed_cases_with_verified_hashes": completed_cases,
               "planned_method_problems": config["planned_method_problems"], "planned_refinements": config["planned_refinements"],
               "preserved_original_fits": config["original_fit_count"], "refinement_states": dict(counters),
               "censored_or_not_recorded": incomplete, "aggregates": aggregates,
               "criteria_failure_counts": dict(Counter(key for row in comparisons for key, passed in row["criteria"].items() if not passed)),
               "counterexamples": [{"case": row["case"], "method": row["method"], "failed_criteria": [key for key, passed in row["criteria"].items() if not passed]}
                                   for row in comparisons if not row["all_predeclared_criteria_met"]],
               "comparisons": comparisons, "production_default_changed": False, "original_outputs_replaced": False,
               "scope": config["scope"], "stationarity_or_global_optimality_certified": False}
    with tempfile.NamedTemporaryFile(mode="w", dir=HERE, delete=False) as handle:
        handle.write(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    os.replace(handle.name, HERE / "summary.json")
    print(json.dumps({key: value for key, value in summary.items() if key not in {"comparisons", "censored_or_not_recorded"}}, indent=2))


if __name__ == "__main__":
    main()
