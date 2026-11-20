"""Refine every frozen coarse endpoint once without discarding original candidates."""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import inspect
import json
import os
from pathlib import Path
import tempfile
import time

import numpy as np
from scipy.optimize import fmin_l_bfgs_b

from numerical_model import make_model, neighborhood, validate_expected_spectrum


HERE = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def atomic_json(path, value):
    with tempfile.NamedTemporaryFile(mode="w", dir=path.parent, delete=False) as handle:
        handle.write(json.dumps(value, indent=2, allow_nan=False) + "\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(handle.name, path)


def run(index):
    manifest_path = HERE / "manifest.json"
    config = json.loads(manifest_path.read_text())
    case = config["cases"][index]
    for name, expected in config["source_sha256"].items():
        if digest(HERE / name) != expected:
            raise ValueError("Frozen diagnostic source changed")
    guard = Path(inspect.getsourcefile(validate_expected_spectrum))
    if digest(guard) != config["likelihood_domain_guard_sha256"]:
        raise ValueError("Likelihood guard changed")
    input_path = HERE / case["input"]
    if digest(input_path) != case["input_sha256"]:
        raise ValueError("Input spectrum changed")
    originals = HERE / "original_results" / case["name"]
    for name, expected in case["original_output_sha256"].items():
        if digest(originals / name) != expected:
            raise ValueError("Original endpoint record changed")
    with np.load(input_path, allow_pickle=False) as data:
        spectrum = data["spectrum"]
    output = HERE / "results" / case["name"]
    output.mkdir(parents=True, exist_ok=False)
    identity = {"case_name": case["name"], "manifest_sha256": digest(manifest_path),
                "parent_manifest_sha256": config["parent_manifest_sha256"], "input_sha256": case["input_sha256"]}
    atomic_json(output / "environment.json", {**identity, "versions": {
        name: importlib.metadata.version(name) for name in ("moments-popgen", "dadi", "scipy", "numpy")},
        "likelihood_guard_sha256": digest(guard)})
    lower, upper = np.log(np.asarray(config["bounds"]))
    bounds = list(zip(lower, upper))
    for method in case["methods"]:
        objective = make_model(method, spectrum, config)
        for restart, start in enumerate(case["refinement_starts"][method]):
            path = output / f"{method}_refined_{restart:02d}.json"
            record = {**identity, "method": method, "arm": "refined", "restart": restart,
                      "start_log_parameters": start["log_parameters"],
                      "original_random_start_log_parameters": case["log_uniform_starts"][restart],
                      "coarse_endpoint_sha256": start["coarse_endpoint_sha256"], "state": "RUNNING"}
            atomic_json(path, record)
            began = time.monotonic()
            try:
                initial_objective = objective(np.asarray(start["log_parameters"]))
                original = json.loads((originals / f"{method}_coarse_{restart:02d}.json").read_text())
                if not np.isclose(initial_objective, original["objective"], rtol=0, atol=1e-8):
                    raise ValueError("Coarse endpoint objective does not reproduce")
                point, _, optimizer = fmin_l_bfgs_b(
                    objective, np.asarray(start["log_parameters"]), bounds=bounds, approx_grad=True,
                    epsilon=config["epsilon"], **config["controls"][method], **config["shared_controls"])
                record.update(state="COMPLETED", parameters=np.exp(point).tolist(),
                              initial_objective=initial_objective, objective=objective(point),
                              optimizer={"converged": bool(optimizer["warnflag"] == 0), "status": int(optimizer["warnflag"]),
                                         "message": str(optimizer["task"]), "evaluations": int(optimizer["funcalls"]),
                                         "iterations": int(optimizer["nit"]), "gradient": np.asarray(optimizer["grad"]).tolist()},
                              neighborhood=neighborhood(objective, point, config))
            except Exception as error:
                record.update(state="FAILED", failure={"type": type(error).__name__, "message": str(error)})
            record["seconds"] = time.monotonic() - began
            atomic_json(path, record)
            print(json.dumps({key: record.get(key) for key in ("case_name", "method", "restart", "state", "objective", "seconds")}), flush=True)
    if digest(input_path) != case["input_sha256"] or digest(manifest_path) != identity["manifest_sha256"]:
        raise ValueError("Diagnostic inputs changed during execution")
    for name, expected in config["source_sha256"].items():
        if digest(HERE / name) != expected:
            raise ValueError("Source changed during execution")
    for name, expected in case["original_output_sha256"].items():
        if digest(originals / name) != expected:
            raise ValueError("Original endpoint changed during execution")
    if digest(guard) != config["likelihood_domain_guard_sha256"]:
        raise ValueError("Likelihood guard changed during execution")
    atomic_json(output / "completion.json", {**identity, "planned_refinements": case["planned_refinements"],
                "output_sha256": {path.name: digest(path) for path in sorted(output.glob("*.json"))},
                "input_and_source_hashes_verified": True})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-index", required=True, type=int)
    run(parser.parse_args().case_index)
