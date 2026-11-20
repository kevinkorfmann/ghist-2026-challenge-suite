"""Execute a frozen paired coarse/fine optimizer case with durable fit records."""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import inspect
from itertools import product
import json
from pathlib import Path
import time

import dadi
import moments
import numpy as np
from scipy.optimize import fmin_l_bfgs_b

from ghist_inference.checkpoint import atomic_json
from ghist_inference.validation import validate_expected_spectrum


HERE = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def make_model(method, spectrum, config):
    package = moments if method == "moments" else dadi
    data = package.Spectrum(spectrum).fold()
    grids = [int(data.sample_sizes[0]) + offset for offset in config["dadi_grid_offsets"]]
    extrapolated = dadi.Numerics.make_extrap_log_func(dadi.Demographics1D.two_epoch)

    def model(point):
        parameters = np.exp(point)
        expected = moments.Demographics1D.two_epoch(parameters, data.sample_sizes) if method == "moments" else extrapolated(parameters, data.sample_sizes, grids)
        return validate_expected_spectrum(expected.fold(), data)

    def objective(point):
        value = -float(package.Inference.ll_multinom(model(point), data))
        if not np.isfinite(value):
            raise ValueError("Nonfinite fixed-domain objective")
        return value

    return objective


def neighborhood(objective, point, config):
    lower, upper = np.log(np.asarray(config["bounds"]))
    center = objective(point)
    records, slopes = [], {}
    for step in config["finite_neighborhood"]["log_steps"]:
        values = {}
        for direction in product((-1, 0, 1), repeat=2):
            if direction == (0, 0):
                continue
            probe = np.clip(point + step * np.asarray(direction), lower, upper)
            value = objective(probe)
            values[direction] = (probe, value)
            records.append({"log_step": step, "direction": direction,
                            "actual_log_displacement": (probe - point).tolist(), "objective_change": value - center})
        dimensions = []
        for index in range(2):
            minus = tuple(-int(i == index) for i in range(2))
            plus = tuple(int(i == index) for i in range(2))
            lo, f_lo = values[minus]
            hi, f_hi = values[plus]
            h_minus, h_plus = point[index] - lo[index], hi[index] - point[index]
            symmetric = bool(h_minus > 0 and h_plus > 0 and np.isclose(h_minus, h_plus, rtol=1e-8, atol=1e-15))
            dimensions.append({"actual_lower_step": h_minus, "actual_upper_step": h_plus,
                               "symmetric": symmetric,
                               "centered_slope": (f_hi - f_lo) / (h_minus + h_plus) if symmetric else None,
                               "forward_slope": (f_hi - center) / h_plus if h_plus > 0 else None,
                               "backward_slope": (center - f_lo) / h_minus if h_minus > 0 else None})
        slopes[str(step)] = dimensions
    improvement = max(0.0, -min(row["objective_change"] for row in records))
    return {"maximum_tested_improvement": improvement,
            "passes_finite_neighborhood_check": improvement <= config["finite_neighborhood"]["maximum_objective_improvement"],
            "probes": records, "axis_slopes": slopes,
            "scope": "feasible finite-neighborhood probes; no KKT/global stationarity certification"}


def run(index):
    manifest_path = HERE / "manifest.json"
    config = json.loads(manifest_path.read_text())
    case = config["cases"][index]
    for filename, expected in config["source_sha256"].items():
        if digest(HERE / filename) != expected:
            raise ValueError("Frozen diagnostic source changed")
    guard = Path(inspect.getsourcefile(validate_expected_spectrum))
    if digest(guard) != config["likelihood_domain_guard_sha256"]:
        raise ValueError("Likelihood guard differs from frozen protocol")
    path = HERE / case["input"]
    if digest(path) != case["input_sha256"]:
        raise ValueError("Frozen spectrum input changed")
    with np.load(path, allow_pickle=False) as data:
        spectrum = data["spectrum"]
    out = HERE / "results" / case["name"]
    out.mkdir(parents=True, exist_ok=False)
    identity = {"case_name": case["name"], "manifest_sha256": digest(manifest_path), "input_sha256": case["input_sha256"]}
    atomic_json(out / "environment.json", {**identity, "versions": {
        name: importlib.metadata.version(name) for name in ("moments-popgen", "dadi", "scipy", "numpy")},
        "likelihood_guard": str(guard), "likelihood_guard_sha256": digest(guard)})
    lower, upper = np.log(np.asarray(config["bounds"]))
    bounds = list(zip(lower, upper))
    objectives = {method: make_model(method, spectrum, config) for method in case["methods"]}
    for restart, start in enumerate(case["log_uniform_starts"]):
        for method, objective in objectives.items():
            for arm, epsilon in config["arms"].items():
                target = out / f"{method}_{arm}_{restart:02d}.json"
                record = {**identity, "method": method, "arm": arm, "restart": restart,
                          "start_log_parameters": start, "state": "RUNNING"}
                atomic_json(target, record)
                began = time.monotonic()
                try:
                    point, _, optimizer = fmin_l_bfgs_b(
                        objective, np.asarray(start), bounds=bounds, approx_grad=True,
                        epsilon=epsilon, **config["controls"][method], **config["shared_controls"])
                    record.update(state="COMPLETED", parameters=np.exp(point).tolist(), objective=objective(point),
                                  optimizer={"converged": bool(optimizer["warnflag"] == 0), "status": int(optimizer["warnflag"]),
                                             "message": str(optimizer["task"]), "evaluations": int(optimizer["funcalls"]),
                                             "iterations": int(optimizer["nit"]), "gradient": np.asarray(optimizer["grad"]).tolist()},
                                  neighborhood=neighborhood(objective, point, config))
                except Exception as error:
                    record.update(state="FAILED", failure={"type": type(error).__name__, "message": str(error)})
                record["seconds"] = time.monotonic() - began
                atomic_json(target, record)
                print(json.dumps({k: record.get(k) for k in ("case_name", "method", "arm", "restart", "state", "objective", "seconds")}), flush=True)
    if digest(path) != case["input_sha256"] or digest(manifest_path) != identity["manifest_sha256"]:
        raise ValueError("Inputs changed during evaluation")
    for filename, expected in config["source_sha256"].items():
        if digest(HERE / filename) != expected:
            raise ValueError("Source changed during evaluation")
    atomic_json(out / "completion.json", {**identity, "planned_fits": case["planned_fits"],
                "output_sha256": {p.name: digest(p) for p in sorted(out.glob("*.json"))},
                "input_and_source_hashes_verified": True})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-index", required=True, type=int)
    run(parser.parse_args().case_index)
