"""Freeze all available complete spectra without reading truth parameters."""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
SEED = 2026113010


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def spectrum_digest(values):
    values = np.asarray(values, dtype="<f8", order="C")
    return hashlib.sha256(str(values.shape).encode() + values.tobytes()).hexdigest()


def main():
    target = HERE / "manifest.json"
    if target.exists():
        raise FileExistsError("Manifest already frozen; do not expand it during execution")
    paths = sorted((HERE.parent / "results").glob("case_*"))
    eligible = []
    for path in paths:
        if not (path / "spectrum.npz").is_file() or not (path / "moments.json").is_file():
            continue
        identity = json.loads((path / "moments.json").read_text())["identity"]
        methods = ["moments", "dadi"] if identity["scenario"] == "matched" else ["moments"]
        if all((path / f"{method}.json").is_file() for method in methods):
            eligible.append((path, identity, methods))
    inputs = HERE / "inputs"
    inputs.mkdir(exist_ok=False)
    records = []

    def add(name, spectrum, identity, methods, sources, length, dataset_code, case_index):
        spectrum = np.asarray(spectrum, dtype=float)
        if spectrum.shape != (115,) or not np.isfinite(spectrum).all() or np.any(spectrum < 0):
            raise ValueError(f"Invalid source spectrum: {name}")
        path = inputs / f"{name}.npz"
        np.savez_compressed(path, spectrum=spectrum)
        rng = np.random.default_rng(np.random.SeedSequence([SEED, dataset_code, case_index]))
        starts = rng.uniform(np.log([0.02, 1e-4]), np.log([0.999, 5.0]), size=(8, 2))
        records.append({
            "name": name, "methods": methods, "input": str(path.relative_to(HERE)),
            "input_sha256": digest(path), "spectrum_sha256": spectrum_digest(spectrum),
            "source_identity": identity, "source_sha256": sources,
            "sequence_length": length, "mutation_rate": 3.42e-8,
            "seed_sequence_entropy": [SEED, dataset_code, case_index],
            "log_uniform_starts": starts.tolist(), "planned_fits": len(methods) * 16,
        })

    for path, identity, methods in eligible:
        with np.load(path / "spectrum.npz", allow_pickle=False) as data:
            blocks = np.asarray(data["spectrum"], dtype=float)
            metadata = json.loads(str(data["metadata"]))
        if spectrum_digest(blocks) != metadata["spectrum_sha256"]:
            raise ValueError("Original block spectrum digest mismatch")
        if any(metadata.get(key) != value for key, value in identity.items()):
            raise ValueError("Original spectrum/fit identity mismatch")
        for method in methods:
            if json.loads((path / f"{method}.json").read_text())["identity"] != identity:
                raise ValueError("Original fit identities differ")
        files = [path / "spectrum.npz", path.parent / "protocol.json", *(path / f"{m}.json" for m in methods)]
        add(f"100mb_{path.name}", blocks.sum(axis=0), identity, methods,
            {str(p.relative_to(ROOT)): digest(p) for p in files}, 100_000_000, 100, identity["case"])
    recovery_path = ROOT / "experiments/scientific_2026_11_20/bottleneck/recovery.json"
    recovery = json.loads(recovery_path.read_text())
    if recovery["complete"] is not True or len(recovery["cases"]) != 24:
        raise ValueError("Require the entire original 24-case 10 Mb panel")
    identities = set()
    for index, case in enumerate(recovery["cases"]):
        identity = {key: case[key] for key in ("point", "replicate", "scenario", "seed")}
        key = (identity["point"], identity["replicate"], identity["scenario"])
        if key in identities:
            raise ValueError("Duplicate original 10 Mb case")
        identities.add(key)
        methods = ["moments", "dadi"] if identity["scenario"] == "matched" else ["moments"]
        add(f"10mb_case_{index:02d}", case["spectrum"], identity, methods,
            {str(recovery_path.relative_to(ROOT)): digest(recovery_path)}, 10_000_000, 10, index)
    manifest = {
        "schema": "bottleneck_paired_epsilon_extension_v1",
        "frozen_at_utc": datetime.now(timezone.utc).isoformat(), "cases": records,
        "planned_cases": len(records), "planned_method_problems": sum(len(c["methods"]) for c in records),
        "planned_fits": sum(c["planned_fits"] for c in records),
        "full_length_panel_planned_cases": 24,
        "full_length_cases_available_at_freeze": [identity["case"] for _, identity, _ in eligible],
        "full_length_cases_unavailable_at_freeze": sorted(set(range(24)) - {identity["case"] for _, identity, _ in eligible}),
        "earlier_10mb_cases": 24, "fresh_seed_root": SEED, "restarts": 8,
        "arms": {"coarse": 1e-3, "fine": 1e-5},
        "bounds": [[0.02, 1e-4], [0.999, 5.0]], "folded": True,
        "dadi_grid_offsets": [20, 30, 40],
        "controls": {"moments": {"maxfun": 15000, "maxiter": 500}, "dadi": {"maxfun": 500, "maxiter": 15000}},
        "shared_controls": {"pgtol": 1e-5, "factr": 1e7, "maxls": 20},
        "finite_neighborhood": {"log_steps": [1e-3, 1e-4, 1e-5], "directions": "all eight nonzero {-1,0,1}^2 axis/diagonal directions, clipped to feasible bounds", "maximum_objective_improvement": 1e-4, "scope": "finite-neighborhood check only; not KKT or stationarity certification"},
        "comparison_criteria": {"maximum_best_objective_deterioration": 1e-3, "no_loss_of_finite_fits": True, "selected_best_reports_converged": True, "minimum_restarts_within_0_1_of_best": 2},
        "source_sha256": {name: digest(HERE / name) for name in ("prepare.py", "run_case.py", "run.sbatch")},
        "likelihood_domain_guard_sha256": digest(ROOT / "src/ghist_inference/validation.py"),
        "allocation": {"cpu_per_case": 1, "maximum_concurrent_cases": 4, "wallclock_seconds_per_case": 600},
        "truth_used_for_fitting_or_selection": False, "original_outputs_replaced": False,
        "limitations": ["Availability of full-length results can favor easier configurations", "Diagnostic parameter panel, not random prior coverage", "Same-spectrum optimizer comparison, not new biological replication", "No source default promotion from a convergence flag alone"],
    }
    target.write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")
    print(json.dumps({key: manifest[key] for key in ("planned_cases", "planned_method_problems", "planned_fits", "full_length_cases_available_at_freeze", "full_length_cases_unavailable_at_freeze")}))


if __name__ == "__main__":
    main()
