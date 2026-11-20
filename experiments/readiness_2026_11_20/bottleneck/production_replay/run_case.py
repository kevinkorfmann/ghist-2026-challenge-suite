"""Run frozen production fitters and compare all retained candidates with diagnostics."""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import time

import numpy as np

from ghist_bottleneck.sfs_baseline import fit_two_epoch_spectrum
from ghist_bottleneck.dadi_reference import fit_two_epoch_dadi_spectrum
from ghist_bottleneck.optimization import optimizer_contract, validate_optimizer_contract
from ghist_inference.checkpoint import atomic_json


HERE = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compare(actual, expected, tolerance):
    mismatches = []
    mapping = {(row["arm"], row["restart"]): row for row in actual}
    for reference in expected:
        key = (reference["arm"], reference["restart"])
        record = mapping.get(key)
        reasons = []
        if record is None or record["state"] != "COMPLETED":
            reasons.append("missing or nonfinite endpoint")
        else:
            if abs(record["objective"] - reference["objective"]) > tolerance["objective_absolute"]:
                reasons.append("objective")
            if not np.allclose(record["parameters"], reference["parameters"], rtol=0, atol=tolerance["parameters_absolute"]):
                reasons.append("parameters")
            for production, diagnostic in (("warnflag", "status"), ("evaluations", "evaluations"), ("iterations", "iterations")):
                if record["optimizer"][production] != reference["optimizer"][diagnostic]:
                    reasons.append(production)
        if reasons:
            mismatches.append({"arm": key[0], "restart": key[1], "reasons": reasons})
    return {"expected_candidates": len(expected), "observed_candidates": len(actual),
            "passed": not mismatches and len(actual) == len(expected), "mismatches": mismatches}


def run(index):
    manifest_path = HERE / "manifest.json"
    config = json.loads(manifest_path.read_text())
    for name, expected in config["source_sha256"].items():
        if digest(HERE / name) != expected:
            raise ValueError("Frozen production replay source changed")
    versions = {name: importlib.metadata.version(name) for name in config["runtime_versions"]}
    if versions != config["runtime_versions"]:
        raise ValueError("Replay must use original diagnostic dependency versions")
    case = config["cases"][index]
    input_path = HERE / case["input"]
    if digest(input_path) != case["input_sha256"]:
        raise ValueError("Frozen spectrum changed")
    with np.load(input_path, allow_pickle=False) as archive:
        spectrum = archive["spectrum"]
    out = HERE / "results" / case["name"]
    out.mkdir(parents=True, exist_ok=False)
    identity = {"case": case["name"], "manifest_sha256": digest(manifest_path), "input_sha256": case["input_sha256"]}
    atomic_json(out / "environment.json", {**identity, "versions": versions})
    for engine in case["methods"]:
        target = out / f"{engine}.json"
        record = {**identity, "engine": engine, "state": "RUNNING"}
        atomic_json(target, record)
        started = time.monotonic()
        try:
            if optimizer_contract(engine) != config["optimizer_contracts"][engine]:
                raise ValueError("Production source/configuration contract changed")
            fit = fit_two_epoch_spectrum if engine == "moments" else fit_two_epoch_dadi_spectrum
            arguments = {"seed": case["seed"]} if "seed" in case else {"initial_log_starts": case["log_starts"]}
            result = fit(spectrum, sequence_length=case["sequence_length"], mutation_rate=config["mutation_rate"],
                         restarts=config["restarts"], folded=config["folded"], **arguments)
            validate_optimizer_contract(result, engine=engine, expected=config["optimizer_contracts"][engine], require_converged=False)
            for attempt in result["optimizer_attempts"]:
                if attempt["arm"] in {"coarse", "fine"} and attempt["start_log_parameters"] != case["log_starts"][attempt["restart"]]:
                    raise ValueError("Production initialization differs from locked starts")
            record.update(state="COMPLETED", fit=result)
            if "expected_attempts" in case:
                record["diagnostic_compatibility"] = compare(result["optimizer_attempts"], case["expected_attempts"][engine], config["comparison_tolerance"])
            try:
                validate_optimizer_contract(result, require_agreement=True, require_neighborhood=True)
                record["numerical_gates"] = {"passed": True}
            except ValueError as error:
                record["numerical_gates"] = {"passed": False, "reason": str(error)}
        except Exception as error:
            record.update(state="FAILED", failure=f"{type(error).__name__}: {error}")
        record["seconds"] = time.monotonic() - started
        atomic_json(target, record)
        print(json.dumps({key: record.get(key) for key in ("case", "engine", "state", "seconds", "diagnostic_compatibility", "numerical_gates", "failure")}), flush=True)
    for name, expected in config["source_sha256"].items():
        if digest(HERE / name) != expected:
            raise ValueError("Frozen replay source changed during execution")
    if digest(input_path) != case["input_sha256"] or digest(manifest_path) != identity["manifest_sha256"]:
        raise ValueError("Replay inputs changed during execution")
    atomic_json(out / "completion.json", {**identity, "output_sha256": {path.name: digest(path) for path in out.glob("*.json")}})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-index", type=int, required=True)
    run(parser.parse_args().case_index)
