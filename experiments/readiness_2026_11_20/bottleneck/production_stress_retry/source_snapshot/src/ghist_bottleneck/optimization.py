"""Versioned bounded optimization; candidate status never changes likelihood selection."""
from __future__ import annotations

import hashlib
import importlib.metadata
from itertools import product
from pathlib import Path

import numpy as np
from scipy.optimize import fmin_l_bfgs_b


DEFAULT_STRATEGY = "paired_coarse_fine_refine_v1"
LEGACY_STRATEGY = "legacy_coarse_v1"
STRATEGIES = (DEFAULT_STRATEGY, LEGACY_STRATEGY)
LOWER = np.array([0.02, 1e-4])
UPPER = np.array([0.999, 5.0])


def optimizer_contract(engine, strategy=DEFAULT_STRATEGY):
    if engine not in {"moments", "dadi"} or strategy not in STRATEGIES:
        raise ValueError("Unknown bottleneck engine or optimizer strategy")
    root = Path(__file__).resolve().parents[1]
    sources = [Path(__file__), root / "ghist_bottleneck" / ("sfs_baseline.py" if engine == "moments" else "dadi_reference.py"),
               root / "ghist_inference/validation.py"]
    return {
        "schema": "ghist_bottleneck_optimizer_v1", "engine": engine, "strategy": strategy,
        "optimizer": "scipy.optimize.fmin_l_bfgs_b", "parameterization": "log",
        "bounds": [LOWER.tolist(), UPPER.tolist()],
        "epsilons": {"coarse": 1e-3} if strategy == LEGACY_STRATEGY else {"coarse": 1e-3, "fine": 1e-5, "refined": 1e-5},
        "controls": {"maxfun": 15000 if engine == "moments" else 500,
                     "maxiter": 500 if engine == "moments" else 15000,
                     "pgtol": 1e-5, "factr": 1e7, "maxls": 20},
        "selection": "minimum finite objective regardless status; original candidates precede refinements on exact ties",
        "restart_ancestry": "one origin per original random start; best finite endpoint per origin",
        "neighborhood": {"log_steps": [1e-3, 1e-4, 1e-5], "directions": "eight axis/diagonal directions clipped to bounds",
                         "maximum_objective_improvement": 1e-4, "scope": "finite-neighborhood probes only; no stationarity or global-optimality certification"},
        "source_sha256": {str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest() for path in sources},
        "versions": {name: importlib.metadata.version(name) for name in ("moments-popgen" if engine == "moments" else "dadi", "numpy", "scipy")},
    }


def validate_optimizer_contract(fit, *, engine=None, expected=None, require_converged=True,
                                require_agreement=False, require_neighborhood=False, allow_legacy=False):
    actual = fit.get("optimizer_contract")
    if not isinstance(actual, dict) or fit.get("engine") != actual.get("engine"):
        raise ValueError("Fit lacks verified engine/optimizer identity; regenerate it with the current fitter")
    engine = actual["engine"] if engine is None else engine
    strategy = actual.get("strategy")
    if strategy == LEGACY_STRATEGY and not allow_legacy:
        raise ValueError("Legacy optimizer fit is development-only; regenerate with the current strategy")
    current = optimizer_contract(engine, strategy if allow_legacy else DEFAULT_STRATEGY)
    if actual != current or (expected is not None and expected != current) or actual.get("engine") != engine:
        raise ValueError("Fit optimizer/source/configuration identity differs; regenerate compatible fits and calibration")
    validate_retained_outcomes(fit)
    optimizer = fit.get("optimizer", {})
    converged = optimizer.get("converged")
    if not isinstance(converged, bool) or fit.get("best_restart_converged") is not converged or converged != (optimizer.get("warnflag") == 0):
        raise ValueError("Selected optimizer convergence metadata is missing or inconsistent")
    if require_converged and not converged:
        raise ValueError("Selected optimizer endpoint did not converge")
    if require_agreement and fit.get("optimizer_diagnostics", {}).get("restarts_within_0_1_log_likelihood", 0) < 2:
        raise ValueError("Fewer than two distinct random-start origins agree near the best fit")
    if require_neighborhood and fit.get("optimizer_neighborhood", {}).get("passes_finite_neighborhood_check") is not True:
        raise ValueError("Selected endpoint lacks a passing finite-neighborhood diagnostic")
    return actual


def validate_retained_outcomes(fit):
    """Recompute selection and unique-origin agreement from retained endpoints."""
    attempts = fit.get("optimizer_attempts")
    restarts = fit.get("restarts")
    if not isinstance(attempts, list) or type(restarts) is not int or restarts < 1:
        raise ValueError("Fit lacks retained optimizer outcomes; regenerate it")
    arms = set(fit["optimizer_contract"]["epsilons"])
    identities = [(row.get("restart"), row.get("arm")) for row in attempts]
    if len(set(identities)) != len(identities) or set(identities) != {(origin, arm) for origin in range(restarts) for arm in arms}:
        raise ValueError("Retained optimizer outcomes have duplicate or missing ancestry/stages")
    if any(row.get("state") not in {"COMPLETED", "FAILED", "SKIPPED"} for row in attempts):
        raise ValueError("Optimizer outcomes are incomplete")
    finite = [row for row in attempts if row["state"] == "COMPLETED"]
    if not finite:
        raise ValueError("Fit has no retained finite endpoint")
    for row in finite:
        parameters = np.asarray(row.get("parameters"), dtype=float)
        optimizer = row.get("optimizer", {})
        if (parameters.shape != (2,) or not np.isfinite(parameters).all() or not np.isfinite(row.get("objective", np.nan))
                or type(optimizer.get("converged")) is not bool or optimizer["converged"] != (optimizer.get("warnflag") == 0)):
            raise ValueError("Invalid retained finite endpoint or optimizer status")
    finite.sort(key=lambda row: (row["arm"] == "refined", row["restart"], row["arm"] == "fine"))
    best = min(finite, key=lambda row: row["objective"])
    origin_best = [min((row for row in finite if row["restart"] == origin), key=lambda row: row["objective"])
                   for origin in range(restarts) if any(row["restart"] == origin for row in finite)]
    if (fit.get("best_restart") != best["restart"] or fit.get("best_optimizer_arm") != best["arm"]
            or fit.get("multinomial_log_likelihood") != -best["objective"]
            or fit.get("post_decline_fraction") != best["parameters"][0]
            or fit.get("scaled_time_2N") != best["parameters"][1]
            or fit.get("optimizer") != best["optimizer"] or fit.get("optimizer_neighborhood") != best.get("neighborhood")):
        raise ValueError("Selected optimizer fit differs from the best retained finite endpoint")
    agreement = sum(row["objective"] - best["objective"] <= 0.1 for row in origin_best)
    if (fit.get("successful_restarts") != len(origin_best)
            or fit.get("converged_restarts") != sum(row["optimizer"]["converged"] for row in origin_best)
            or fit.get("optimizer_diagnostics", {}).get("restarts_within_0_1_log_likelihood") != agreement):
        raise ValueError("Optimizer counts disagree with distinct retained random-start origins")
    probe = best.get("neighborhood", {})
    if probe.get("status") == "complete":
        records = probe.get("probes", [])
        expected_probes = {(step, direction) for step in (1e-3, 1e-4, 1e-5)
                           for direction in product((-1, 0, 1), repeat=2) if direction != (0, 0)}
        observed_probes = [(row.get("log_step"), tuple(row.get("direction", []))) for row in records]
        changes = [row.get("objective_change", np.nan) for row in records]
        if len(observed_probes) != 24 or set(observed_probes) != expected_probes or not np.isfinite(changes).all():
            raise ValueError("Selected endpoint has incomplete finite-neighborhood evidence")
        point = np.log(best["parameters"])
        for row in records:
            displacement = np.asarray(row.get("actual_log_displacement"), dtype=float)
            expected_displacement = np.clip(point + row["log_step"] * np.asarray(row["direction"]),
                                            np.log(LOWER), np.log(UPPER)) - point
            if displacement.shape != (2,) or not np.allclose(displacement, expected_displacement, rtol=0, atol=2e-15):
                raise ValueError("Selected neighborhood displacements differ from its endpoint and bounds")
        improvement = max(0.0, -min(changes))
        if (probe.get("maximum_tested_improvement") != improvement
                or probe.get("passes_finite_neighborhood_check") != (improvement <= 1e-4)):
            raise ValueError("Selected neighborhood flag differs from retained raw probes")
    elif probe.get("status") != "failed" or probe.get("passes_finite_neighborhood_check") is not False:
        raise ValueError("Selected endpoint has unknown neighborhood diagnostic status")
    try:
        theta, rate = float(fit["theta"]), float(fit["mutation_rate"])
        exposure = float(fit.get("callable_sequence_length", fit["sequence_length"]))
        size, generations = float(fit["ancestral_effective_size"]), float(fit["generations"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError("Fit lacks physical parameters linked to its selected endpoint") from error
    if (best.get("theta") != theta or not np.isfinite([theta, rate, exposure, size, generations]).all()
            or min(theta, rate, exposure, size, generations) <= 0
            or not np.isclose(size, theta / (4 * rate * exposure), rtol=1e-12, atol=0)
            or not np.isclose(generations, 2 * size * best["parameters"][1], rtol=1e-12, atol=0)):
        raise ValueError("Physical parameters differ from the selected optimizer endpoint and exposure")


def log_starts(restarts, seed, *, initial=None, starts=None):
    """Preserve historical RNG consumption, including an optional first warm start."""
    lower, upper = np.log(LOWER), np.log(UPPER)
    if starts is not None:
        if initial is not None:
            raise ValueError("Explicit log starts and initial estimate are mutually exclusive")
        values = np.asarray(starts, dtype=float)
        if values.shape != (restarts, 2) or not np.isfinite(values).all() or np.any(values < lower) or np.any(values > upper):
            raise ValueError("Explicit log starts must match restarts and lie within bounds")
        return values.copy()
    if initial is not None:
        initial = np.asarray(initial, dtype=float)
        if initial.shape != (2,) or not np.isfinite(initial).all() or np.any(initial < LOWER) or np.any(initial > UPPER):
            raise ValueError("Initial estimate must lie within parameter bounds")
    rng = np.random.default_rng(seed)
    return np.asarray([np.log(initial) if index == 0 and initial is not None
                       else np.log(np.exp(rng.uniform(lower, upper))) for index in range(restarts)])


def neighborhood(objective, point):
    lower, upper = np.log(LOWER), np.log(UPPER)
    center = objective(point)
    records, slopes = [], {}
    for step in (1e-3, 1e-4, 1e-5):
        values = {}
        for direction in product((-1, 0, 1), repeat=2):
            if direction == (0, 0):
                continue
            probe = np.clip(point + step * np.asarray(direction), lower, upper)
            value = objective(probe)
            if not np.isfinite(value):
                raise ValueError("Nonfinite finite-neighborhood objective")
            values[direction] = (probe, value)
            records.append({"log_step": step, "direction": list(direction), "actual_log_displacement": (probe - point).tolist(),
                            "objective_change": float(value - center)})
        dimensions = []
        for index in range(2):
            lo, f_lo = values[tuple(-int(i == index) for i in range(2))]
            hi, f_hi = values[tuple(int(i == index) for i in range(2))]
            h_minus, h_plus = point[index] - lo[index], hi[index] - point[index]
            symmetric = bool(h_minus > 0 and h_plus > 0 and np.isclose(h_minus, h_plus, rtol=1e-8, atol=1e-15))
            dimensions.append({"actual_lower_step": float(h_minus), "actual_upper_step": float(h_plus), "symmetric": symmetric,
                               "centered_slope": float((f_hi - f_lo) / (h_minus + h_plus)) if symmetric else None,
                               "forward_slope": float((f_hi - center) / h_plus) if h_plus > 0 else None,
                               "backward_slope": float((center - f_lo) / h_minus) if h_minus > 0 else None})
        slopes[str(step)] = dimensions
    improvement = max(0.0, -min(row["objective_change"] for row in records))
    return {"status": "complete", "maximum_tested_improvement": improvement,
            "passes_finite_neighborhood_check": improvement <= 1e-4, "probes": records, "axis_slopes": slopes,
            "scope": "feasible finite-neighborhood probes; no KKT/global stationarity certification"}


def run_optimizer(objective, *, engine, starts, strategy=DEFAULT_STRATEGY, evaluate=None, legacy_optimize=None):
    contract = optimizer_contract(engine, strategy)
    starts = np.asarray(starts, dtype=float)
    if starts.ndim != 2 or starts.shape[1] != 2 or len(starts) < 1 or not np.isfinite(starts).all():
        raise ValueError("Optimizer requires finite two-dimensional log starts")
    bounds = list(zip(np.log(LOWER), np.log(UPPER)))
    attempts, candidates = [], []

    def attempt(origin, arm, start):
        record = {"restart": origin, "arm": arm, "start_log_parameters": np.asarray(start).tolist(), "state": "FAILED"}
        attempts.append(record)
        try:
            if strategy == LEGACY_STRATEGY and legacy_optimize is not None:
                point, optimizer = legacy_optimize(start)
            else:
                point, _, optimizer = fmin_l_bfgs_b(objective, np.asarray(start), bounds=bounds,
                                                  approx_grad=True, epsilon=contract["epsilons"][arm], **contract["controls"])
            value = float(objective(point))
            parameters = np.exp(point)
            if (not np.isfinite([*parameters, value]).all()
                    or np.any(np.asarray(point) < np.log(LOWER) - 1e-12)
                    or np.any(np.asarray(point) > np.log(UPPER) + 1e-12)):
                raise ValueError("Nonfinite or out-of-bounds endpoint")
            extra = {} if evaluate is None else evaluate(parameters)
            # Preserve the finite endpoint before optional neighborhood diagnostics.
            record.update(state="COMPLETED", parameters=parameters.tolist(), objective=value, **extra,
                          optimizer={"warnflag": int(optimizer["warnflag"]), "converged": bool(optimizer["warnflag"] == 0),
                                     "message": str(optimizer["task"]), "evaluations": int(optimizer["funcalls"]),
                                     "iterations": int(optimizer.get("nit", 0)), "gradient": np.asarray(optimizer.get("grad", [])).tolist()})
            candidates.append(record)
        except Exception as error:
            record["failure"] = f"{type(error).__name__}: {error}"
            return None
        try:
            record["neighborhood"] = neighborhood(objective, np.asarray(point))
        except Exception as error:
            record["neighborhood"] = {"status": "failed", "passes_finite_neighborhood_check": False,
                                       "failure": f"{type(error).__name__}: {error}"}
        return record

    coarse = []
    for origin, start in enumerate(starts):
        endpoint = attempt(origin, "coarse", start)
        if endpoint is not None:
            coarse.append(endpoint)
        if strategy == DEFAULT_STRATEGY:
            attempt(origin, "fine", start)
    if strategy == DEFAULT_STRATEGY:
        for endpoint in coarse:
            attempt(endpoint["restart"], "refined", np.log(endpoint["parameters"]))
        for origin in sorted(set(range(len(starts))) - {record["restart"] for record in coarse}):
            attempts.append({"restart": origin, "arm": "refined", "state": "SKIPPED", "failure": "No finite coarse endpoint to refine"})
    if not candidates:
        raise RuntimeError(f"All {engine} optimizer attempts failed: {attempts[:3]}")
    best = min(candidates, key=lambda row: row["objective"])
    origin_best = [min((row for row in candidates if row["restart"] == origin), key=lambda row: row["objective"])
                   for origin in range(len(starts)) if any(row["restart"] == origin for row in candidates)]
    arm_summary = {}
    for arm in contract["epsilons"]:
        records = [row for row in attempts if row["arm"] == arm]
        arm_summary[arm] = {"planned": len(starts), "finite": sum(row["state"] == "COMPLETED" for row in records),
                            "reported_converged": sum(row.get("optimizer", {}).get("converged", False) for row in records),
                            "failed": sum(row["state"] == "FAILED" for row in records),
                            "skipped": sum(row["state"] == "SKIPPED" for row in records),
                            "neighborhood_failed": sum(row.get("neighborhood", {}).get("status") == "failed" for row in records)}
    return {"contract": contract, "best": best, "origin_best": origin_best, "attempts": attempts,
            "arm_summary": arm_summary,
            "failed_origins": [{"restart": origin, "reason": "All candidate stages failed"} for origin in range(len(starts))
                               if not any(row["restart"] == origin for row in candidates)]}
