"""Audit the frozen production replay, retaining every planned problem in denominators."""
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

from ghist_bottleneck.optimization import validate_retained_outcomes
from ghist_inference.checkpoint import atomic_json


HERE = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def summarize():
    manifest_path = HERE / "manifest.json"
    config = json.loads(manifest_path.read_text())
    manifest_hash = digest(manifest_path)
    for name, expected in config["source_sha256"].items():
        if digest(HERE / name) != expected:
            raise ValueError(f"Frozen source hash mismatch: {name}")
    loader = importlib.util.spec_from_file_location("frozen_replay", HERE / "run_case.py")
    runner = importlib.util.module_from_spec(loader)
    loader.loader.exec_module(runner)
    rows = []
    completed_cases = 0
    for case in config["cases"]:
        if digest(HERE / case["input"]) != case["input_sha256"]:
            raise ValueError("Frozen input hash mismatch")
        directory = HERE / "results" / case["name"]
        completion_path = directory / "completion.json"
        completion = json.loads(completion_path.read_text()) if completion_path.exists() else None
        if completion:
            if completion["manifest_sha256"] != manifest_hash or completion["input_sha256"] != case["input_sha256"]:
                raise ValueError("Completion input identity mismatch")
            for name, expected in completion["output_sha256"].items():
                if digest(directory / name) != expected:
                    raise ValueError("Output hash mismatch")
            completed_cases += 1
        for engine in case["methods"]:
            path = directory / f"{engine}.json"
            row = {"case": case["name"], "engine": engine,
                   "scope": case.get("scope", "compatibility" if "expected_attempts" in case else "separate_stress"),
                   "planned_candidates": 3 * config["restarts"], "state": "MISSING"}
            if path.exists():
                record = json.loads(path.read_text())
                row.update(state=record["state"], output_sha256=digest(path), seconds=record.get("seconds"), method_attempted=True)
                retained = record.get("fit", {}).get("optimizer_attempts")
                row["retained_attempt_records"] = len(retained) if isinstance(retained, list) else None
                row["attempted_stages_from_retained_records"] = (sum(item["state"] != "SKIPPED" for item in retained)
                                                                 if isinstance(retained, list) else None)
                if record["manifest_sha256"] != manifest_hash or record["input_sha256"] != case["input_sha256"]:
                    raise ValueError("Fit input identity mismatch")
                if record["state"] == "COMPLETED":
                    if completion is None or path.name not in completion["output_sha256"]:
                        row["state"] = "UNATTESTED"
                    else:
                        fit = record["fit"]
                        if fit["optimizer_contract"] != config["optimizer_contracts"][engine]:
                            raise ValueError("Fit engine/strategy/source/configuration mismatch")
                        validate_retained_outcomes(fit)
                        attempts = fit["optimizer_attempts"]
                        selected = next(item for item in attempts if item["restart"] == fit["best_restart"] and item["arm"] == fit["best_optimizer_arm"])
                        candidates = [item for item in attempts if item["state"] == "COMPLETED"]
                        origins = [min((item for item in candidates if item["restart"] == origin), key=lambda item: item["objective"])
                                   for origin in range(config["restarts"]) if any(item["restart"] == origin for item in candidates)]
                        agreement = sum(item["objective"] - selected["objective"] <= .1 for item in origins)
                        row.update(finite_candidates=len(candidates), failed_candidates=sum(item["state"] == "FAILED" for item in attempts),
                                   skipped_candidates=sum(item["state"] == "SKIPPED" for item in attempts),
                                   objective=selected["objective"], selected_arm=selected["arm"], selected_origin=selected["restart"],
                                   selected_converged=selected["optimizer"]["converged"], selected_termination=selected["optimizer"]["message"],
                                   distinct_near_best_origins=agreement, successful_origins=len(origins),
                                   finite_neighborhood_pass=selected["neighborhood"]["passes_finite_neighborhood_check"],
                                   maximum_tested_improvement=selected["neighborhood"].get("maximum_tested_improvement"),
                                   numerical_gates_pass=bool(selected["optimizer"]["converged"] and agreement >= 2
                                                             and selected["neighborhood"]["passes_finite_neighborhood_check"]))
                        coarse = [item["objective"] for item in candidates if item["arm"] == "coarse"]
                        row["selected_minus_best_coarse_objective"] = selected["objective"] - min(coarse) if coarse else None
                        if "expected_attempts" in case:
                            comparison = runner.compare(attempts, case["expected_attempts"][engine], config["comparison_tolerance"])
                            if comparison != record["diagnostic_compatibility"]:
                                raise ValueError("Stored compatibility metrics do not recompute")
                            row["diagnostic_compatibility"] = comparison
                elif "failure" in record:
                    row["failure"] = record["failure"]
            rows.append(row)
    complete = [row for row in rows if row["state"] == "COMPLETED"]
    comparisons = [row for row in complete if "diagnostic_compatibility" in row]
    compatibility = [row for row in rows if row["scope"] == "compatibility"]
    result = {"schema": "bottleneck_production_replay_audit_v1", "manifest_sha256": manifest_hash,
              "audit_source_sha256": digest(Path(__file__)),
              "audit_guard_source_sha256": digest(Path(validate_retained_outcomes.__code__.co_filename)), "planned_cases": config["planned_cases"],
              "completed_cases": completed_cases, "planned_method_problems": config["planned_method_problems"],
              "completed_method_problems": len(complete),
              "attempted_method_problems": sum(row.get("method_attempted", False) for row in rows),
              "planned_optimizer_attempts": config["planned_optimizer_attempts"],
              "retained_optimizer_attempt_records": sum(row.get("retained_attempt_records") or 0 for row in rows),
              "attempted_stages_from_retained_records": sum(row.get("attempted_stages_from_retained_records") or 0 for row in rows),
              "problems_without_retained_stage_counts": [{"case": row["case"], "engine": row["engine"], "state": row["state"]}
                                                        for row in rows if row.get("retained_attempt_records") is None],
              "finite_optimizer_attempts": sum(row["finite_candidates"] for row in complete),
              "endpoint_comparisons_planned": sum(len(expected) for case in config["cases"]
                                                   for expected in case.get("expected_attempts", {}).values()),
              "endpoint_comparisons_matched": sum(row["diagnostic_compatibility"]["expected_candidates"]
                                                  - len(row["diagnostic_compatibility"]["mismatches"]) for row in comparisons),
              "compatibility_planned_method_problems": len(compatibility),
              "compatibility_passed_method_problems": sum(row.get("diagnostic_compatibility", {}).get("passed", False) for row in compatibility),
              "selected_converged_problems": sum(row["selected_converged"] for row in complete),
              "finite_neighborhood_passed_problems": sum(row["finite_neighborhood_pass"] for row in complete),
              "distinct_origin_agreement_passed_problems": sum(row["distinct_near_best_origins"] >= 2 for row in complete),
              "numerical_gates_passed_problems": sum(row["numerical_gates_pass"] for row in complete), "problems": rows,
              "limitations": ["Compatibility and repair development using existing spectra, not independent biological confirmation.",
                              "The separately reported stress spectrum was outside the original 53 numerical problems.",
                              "Finite-neighborhood probes and objective-reduction termination do not certify stationarity or biological correctness.",
                              "Counts conditional on completed problems if any planned problem remains missing, failed or unattested."]}
    atomic_json(HERE / "summary.json", result)
    print(json.dumps({key: value for key, value in result.items() if key not in {"problems", "limitations"}}, indent=2))


if __name__ == "__main__":
    summarize()
