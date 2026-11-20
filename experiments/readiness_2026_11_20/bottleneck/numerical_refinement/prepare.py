"""Freeze one refinement from every coarse endpoint, preserving random-start ancestry."""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np


HERE = Path(__file__).resolve().parent
PARENT = HERE.parent / "numerical_extension"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    target = HERE / "manifest.json"
    if target.exists():
        raise FileExistsError("Refinement already frozen")
    parent_hash = digest(PARENT / "manifest.json")
    if parent_hash != "8c2e403832193cb437e472baecc371aa9f4e6809122df5a027b0d136c16ae7ba":
        raise ValueError("Unexpected parent experiment")
    parent = json.loads((PARENT / "manifest.json").read_text())
    if digest(HERE / "numerical_model.py") != parent["source_sha256"]["run_case.py"]:
        raise ValueError("Model/probe source changed from paired experiment")
    (HERE / "inputs").mkdir(exist_ok=False)
    (HERE / "original_results").mkdir(exist_ok=False)
    cases = []
    for case in parent["cases"]:
        source = PARENT / "results" / case["name"]
        completion = json.loads((source / "completion.json").read_text())
        if completion["manifest_sha256"] != parent_hash or completion["input_sha256"] != case["input_sha256"]:
            raise ValueError("Parent completion identity mismatch")
        for name, expected in completion["output_sha256"].items():
            if digest(source / name) != expected:
                raise ValueError("Parent completed output hash mismatch")
        if digest(PARENT / case["input"]) != case["input_sha256"]:
            raise ValueError("Parent input changed")
        shutil.copyfile(PARENT / case["input"], HERE / case["input"])
        destination = HERE / "original_results" / case["name"]
        shutil.copytree(source, destination)
        starts = {}
        for method in case["methods"]:
            starts[method] = []
            for restart in range(parent["restarts"]):
                path = destination / f"{method}_coarse_{restart:02d}.json"
                coarse = json.loads(path.read_text())
                if (coarse["state"] != "COMPLETED" or coarse["case_name"] != case["name"]
                        or coarse["method"] != method or coarse["restart"] != restart
                        or coarse["start_log_parameters"] != case["log_uniform_starts"][restart]):
                    raise ValueError("Coarse endpoint missing or ancestry changed")
                parameters = np.asarray(coarse["parameters"], dtype=float)
                if not np.isfinite(parameters).all() or np.any(parameters <= 0):
                    raise ValueError("Coarse parameters invalid")
                starts[method].append({"restart": restart, "log_parameters": np.log(parameters).tolist(),
                                       "coarse_endpoint_sha256": digest(path)})
        cases.append({**case, "planned_refinements": len(case["methods"]) * parent["restarts"],
                      "refinement_starts": starts,
                      "original_output_sha256": {path.name: digest(path) for path in sorted(destination.glob("*.json"))}})
    manifest = {
        "schema": "bottleneck_coarse_endpoint_refinement_v1",
        "frozen_at_utc": datetime.now(timezone.utc).isoformat(),
        "parent_manifest_sha256": parent_hash,
        "cases": cases, "planned_cases": len(cases),
        "planned_method_problems": parent["planned_method_problems"],
        "planned_refinements": sum(case["planned_refinements"] for case in cases),
        "original_fit_count": parent["planned_fits"], "restarts": parent["restarts"],
        "epsilon": 1e-5,
        **{key: parent[key] for key in ("bounds", "folded", "dadi_grid_offsets", "controls", "shared_controls", "finite_neighborhood", "likelihood_domain_guard_sha256")},
        "selection": "minimum finite objective across original coarse, original fine and coarse-endpoint refinement; exact ties preserve original candidate; no convergence-status filtering",
        "restart_agreement": "minimum objective over candidates sharing each original random-start index; a refinement and its ancestor count once",
        "criteria": {"all_planned_refinements_complete": True, "no_objective_degradation_vs_original_union": True,
                     "selected_best_reports_converged": True, "minimum_distinct_origins_within_0_1": 2,
                     "selected_best_passes_existing_finite_neighborhood": True},
        "allocation": {"cpu_per_case": 1, "maximum_concurrent_cases": 4, "wallclock_seconds_per_case": 600},
        "source_sha256": {name: digest(HERE / name) for name in ("prepare.py", "run_case.py", "numerical_model.py", "run.sbatch")},
        "runtime_scope": "package versions recorded; external package sources not hash-locked",
        "truth_used_for_fitting_or_selection": False, "original_outputs_replaced": False,
        "production_default_changed": False,
        "scope": "Prespecified repair development on all 53 existing problems; not new biological recovery or independent confirmation; no subsequent parameter search authorized by this protocol",
    }
    target.write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")
    print(json.dumps({key: manifest[key] for key in ("planned_cases", "planned_method_problems", "planned_refinements", "original_fit_count")}))
    print("manifest_sha256", digest(target))


if __name__ == "__main__":
    main()
