"""Isolated optimizer diagnostic; never replaces a locked panel result.

All starts and variants are fixed before execution. No simulation truth is
used for initialization, optimizer comparison or selection.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import inspect
import json
from pathlib import Path
import time

import dadi
import moments
import numpy as np
from scipy.optimize import fmin_l_bfgs_b, minimize

from ghist_inference.checkpoint import atomic_json, load_spectrum_checkpoint
from ghist_inference.validation import validate_expected_spectrum


SEED = 2026113008
RESTARTS = 8
LOWER = np.array([0.02, 1e-4])
UPPER = np.array([0.999, 5.0])
GRADIENT_STEPS = (1e-3, 1e-4, 1e-5)
VARIANTS = (
    {"name": "lbfgsb_epsilon_1e-3", "optimizer": "lbfgsb", "epsilon": 1e-3},
    {"name": "lbfgsb_epsilon_1e-5", "optimizer": "lbfgsb", "epsilon": 1e-5},
    {"name": "powell", "optimizer": "powell", "xtol": 1e-5, "ftol": 1e-10},
)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def make_model(method, spectrum):
    package = moments if method == "moments" else dadi
    data = package.Spectrum(spectrum).fold()
    grids = [int(data.sample_sizes[0]) + offset for offset in (20, 30, 40)]
    extrapolated = dadi.Numerics.make_extrap_log_func(dadi.Demographics1D.two_epoch)

    def model(log_parameters):
        parameters = np.exp(log_parameters)
        if method == "moments":
            expected = moments.Demographics1D.two_epoch(parameters, data.sample_sizes)
        else:
            expected = extrapolated(parameters, data.sample_sizes, grids)
        return validate_expected_spectrum(expected.fold(), data)

    def objective(log_parameters):
        likelihood = float(package.Inference.ll_multinom(model(log_parameters), data))
        if not np.isfinite(likelihood):
            raise ValueError("Nonfinite likelihood")
        return -likelihood

    def theta(log_parameters):
        return float(package.Inference.optimal_sfs_scaling(model(log_parameters), data))

    return objective, theta


def local_differences(objective, point):
    """Fixed finite-difference checks; no follow-up point is optimized or chosen."""
    center = objective(point)
    checks = {}
    for step in GRADIENT_STEPS:
        dimensions = []
        for index in range(2):
            low, high = point.copy(), point.copy()
            low[index] = max(low[index] - step, np.log(LOWER[index]))
            high[index] = min(high[index] + step, np.log(UPPER[index]))
            f_low, f_high = objective(low), objective(high)
            dimensions.append({
                "central_slope": float((f_high - f_low) / (high[index] - low[index])),
                "forward_slope": float((f_high - center) / (high[index] - point[index]))
                if high[index] > point[index] else None,
                "backward_slope": float((center - f_low) / (point[index] - low[index]))
                if point[index] > low[index] else None,
                "lower_objective_change": f_low - center,
                "upper_objective_change": f_high - center,
                "centered": bool(low[index] < point[index] < high[index]),
            })
        checks[str(step)] = dimensions
    return checks


def run(args):
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    case = args.case.resolve()
    blocks, metadata = load_spectrum_checkpoint(case / "spectrum.npz", expected={"case": 8})
    spectrum = blocks.sum(axis=0)
    references = {method: json.loads((case / f"{method}.json").read_text()) for method in ("moments", "dadi")}
    if digest(case / "panel_protocol.json") != metadata["protocol_sha256"]:
        raise ValueError("Original protocol digest differs from the spectrum identity")
    for record in references.values():
        if any(metadata.get(key) != value for key, value in record["identity"].items()):
            raise ValueError("Original fit and spectrum identities differ")
    rng = np.random.default_rng(SEED)
    starts = rng.uniform(np.log(LOWER), np.log(UPPER), size=(RESTARTS, 2))
    source_paths = {
        "diagnostic": Path(__file__).resolve(),
        "moments_demography": Path(inspect.getsourcefile(moments.Demographics1D.two_epoch)),
        "dadi_demography": Path(inspect.getsourcefile(dadi.Demographics1D.two_epoch)),
        "scipy_lbfgsb": Path(inspect.getsourcefile(fmin_l_bfgs_b)),
        "scipy_minimize": Path(inspect.getsourcefile(minimize)),
        "likelihood_domain_guard": Path(inspect.getsourcefile(validate_expected_spectrum)),
        "spectrum_checkpoint_loader": Path(inspect.getsourcefile(load_spectrum_checkpoint)),
    }
    inputs = {str(case / name): digest(case / name) for name in ("spectrum.npz", "moments.json", "dadi.json", "panel_protocol.json")}
    sources = {name: {"path": str(path), "sha256": digest(path)} for name, path in source_paths.items()}
    protocol = {
        "schema": "case08_numerical_diagnostic_v1", "fresh_seed": SEED,
        "restarts_per_variant": RESTARTS, "planned_fits": 48,
        "variants": VARIANTS, "methods": ["moments", "dadi"],
        "log_uniform_starts": starts.tolist(), "parameter_bounds": [LOWER.tolist(), UPPER.tolist()],
        "fixed_model": "folded two-epoch multinomial SFS; same fixed likelihood domain as original panel",
        "dadi_grid_offsets": [20, 30, 40],
        "lbfgsb_controls": {"pgtol": 1e-5, "factr": 1e7, "maxfun": 1000, "maxiter": 500, "maxls": 20},
        "powell_controls": {"maxfev": 1000, "maxiter": 500},
        "local_difference_steps": GRADIENT_STEPS,
        "allocated_cpus": 1, "slurm_wallclock_limit_seconds": 600,
        "input_sha256": inputs, "source_sha256": sources,
        "spectrum_sha256": metadata["spectrum_sha256"],
        "original_identity": references["moments"]["identity"],
        "versions": {name: importlib.metadata.version(name) for name in ("moments-popgen", "dadi", "scipy", "numpy")},
        "selection": "Case selected for abnormal selected optimizers; diagnostic, not an independent recovery replicate",
        "truth_used": False, "panel_results_replaced": False,
        "scope": "Compare numerical termination and objective agreement; no biological adequacy or promotion claim",
    }
    atomic_json(out / "protocol.json", protocol)
    models = {method: make_model(method, spectrum) for method in ("moments", "dadi")}
    original = {}
    for method, (objective, _) in models.items():
        fit = references[method]["fit"]
        point = np.log([fit["post_decline_fraction"], fit["scaled_time_2N"]])
        original[method] = {
            "parameters": np.exp(point).tolist(), "stored_objective": -fit["multinomial_log_likelihood"],
            "recomputed_objective": objective(point), "optimizer": fit["optimizer"],
            "local_differences": local_differences(objective, point),
        }
    atomic_json(out / "original_selected_diagnostics.json", original)
    results = []
    bounds = list(zip(np.log(LOWER), np.log(UPPER)))
    for restart, start in enumerate(starts):
        for method, (objective, theta) in models.items():
            for variant in VARIANTS:
                began = time.monotonic()
                row = {"method": method, "variant": variant["name"], "restart": restart, "start_log_parameters": start.tolist()}
                try:
                    if variant["optimizer"] == "lbfgsb":
                        point, _, optimizer = fmin_l_bfgs_b(
                            objective, start.copy(), bounds=bounds, approx_grad=True,
                            epsilon=variant["epsilon"], pgtol=1e-5, factr=1e7,
                            maxfun=1000, maxiter=500, maxls=20,
                        )
                        status = {"converged": bool(optimizer["warnflag"] == 0),
                                  "status": int(optimizer["warnflag"]), "message": str(optimizer["task"]),
                                  "evaluations": int(optimizer["funcalls"]), "iterations": int(optimizer["nit"]),
                                  "optimizer_gradient": np.asarray(optimizer["grad"]).tolist()}
                    else:
                        fitted = minimize(objective, start.copy(), method="Powell", bounds=bounds,
                                          options={"maxfev": 1000, "maxiter": 500,
                                                   "xtol": variant["xtol"], "ftol": variant["ftol"]})
                        point = fitted.x
                        status = {"converged": bool(fitted.success), "status": int(fitted.status),
                                  "message": str(fitted.message), "evaluations": int(fitted.nfev),
                                  "iterations": int(fitted.nit)}
                    row.update(parameters=np.exp(point).tolist(), objective=objective(point),
                               theta=theta(point), optimizer=status,
                               local_differences=local_differences(objective, point))
                except Exception as error:
                    row["failure"] = {"type": type(error).__name__, "message": str(error)}
                row["seconds"] = time.monotonic() - began
                results.append(row)
                atomic_json(out / "fits.json", results)
                print(json.dumps({k: row.get(k) for k in ("method", "variant", "restart", "objective", "seconds", "failure")}), flush=True)
    for filename, expected in inputs.items():
        if digest(filename) != expected:
            raise ValueError(f"Original input changed: {filename}")
    for source in sources.values():
        if digest(source["path"]) != source["sha256"]:
            raise ValueError(f"Source changed: {source['path']}")
    atomic_json(out / "completion.json", {
        "planned_fits": 48, "recorded_fits": len(results), "input_and_source_hashes_verified": True,
        "protocol_sha256": digest(out / "protocol.json"),
        "output_sha256": {name: digest(out / name) for name in ("fits.json", "original_selected_diagnostics.json")},
        "panel_results_replaced": False,
    })


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    run(parser.parse_args())
