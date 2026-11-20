"""Freeze production replay against existing numerical evidence, without truth access."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np

from ghist_bottleneck.optimization import log_starts, optimizer_contract


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(diagnostics, stress):
    if (HERE / "manifest.json").exists():
        raise FileExistsError("Replay is already frozen")
    extension = diagnostics / "numerical_extension"
    refinement = diagnostics / "numerical_refinement"
    parent = json.loads((extension / "manifest.json").read_text())
    records = {case["name"]: case for case in parent["cases"]}
    cases = []
    (HERE / "inputs").mkdir(exist_ok=False)
    for name in ("100mb_case_08", "100mb_case_02", "10mb_case_16", "10mb_case_20"):
        old = records[name]
        target = HERE / "inputs" / f"{name}.npz"
        shutil.copyfile(extension / old["input"], target)
        expected, source_hashes = {}, {}
        for method in old["methods"]:
            expected[method] = []
            for arm in ("coarse", "fine", "refined"):
                source = extension if arm != "refined" else refinement
                for restart in range(8):
                    path = source / "results" / name / f"{method}_{arm}_{restart:02d}.json"
                    expected[method].append(json.loads(path.read_text()))
                    source_hashes[f"{source.name}/results/{name}/{path.name}"] = digest(path)
        cases.append({"name": name, "input": str(target.relative_to(HERE)), "input_sha256": digest(target),
                      "sequence_length": old["sequence_length"], "methods": old["methods"],
                      "log_starts": old["log_uniform_starts"], "initialization": "identical explicit diagnostic log starts",
                      "expected_attempts": expected, "reference_sha256": source_hashes})
    with np.load(stress, allow_pickle=False) as archive:
        spectrum = archive["spectrum"].sum(axis=0)
    target = HERE / "inputs/100mb_case_16_new_stress.npz"
    np.savez_compressed(target, spectrum=spectrum)
    cases.append({"name": "100mb_case_16_new_stress", "input": str(target.relative_to(HERE)),
                  "input_sha256": digest(target), "sequence_length": 100_000_000,
                  "methods": ["moments", "dadi"], "seed": 2026112951,
                  "log_starts": log_starts(8, 2026112951).tolist(),
                  "initialization": "original seeded random starts; no truth initialization",
                  "reference_sha256": {str(stress): digest(stress)},
                  "scope": "one separately reported new stress spectrum outside the original53 numerical problems"})
    sources = ["ghist_bottleneck/__init__.py", "ghist_bottleneck/optimization.py", "ghist_bottleneck/sfs_baseline.py",
               "ghist_bottleneck/dadi_reference.py", "ghist_inference/__init__.py", "ghist_inference/validation.py",
               "ghist_inference/checkpoint.py", "ghist_inference/vcf.py"]
    source_hashes = {}
    for name in sources:
        target = HERE / "source_snapshot/src" / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / "src" / name, target)
        source_hashes[str(target.relative_to(HERE))] = digest(target)
    versions = {"moments-popgen": "1.6.0", "dadi": "2.4.4", "numpy": "2.4.6", "scipy": "1.16.3"}
    contracts = {engine: optimizer_contract(engine) for engine in ("moments", "dadi")}
    for contract in contracts.values():
        contract["versions"] = {package: versions[package] for package in contract["versions"]}
    manifest = {"schema": "bottleneck_production_replay_v1", "cases": cases, "planned_cases": 5,
                "planned_method_problems": 10, "planned_optimizer_attempts": 240,
                "mutation_rate": 3.42e-8, "restarts": 8, "folded": True,
                "runtime_versions": versions, "optimizer_contracts": contracts,
                "source_sha256": {**source_hashes, **{name: digest(HERE / name) for name in ("prepare.py", "run_case.py", "run.sbatch")}},
                "extension_manifest_sha256": digest(extension / "manifest.json"),
                "refinement_manifest_sha256": digest(refinement / "manifest.json"),
                "comparison_tolerance": {"objective_absolute": 1e-8, "parameters_absolute": 1e-10,
                                         "optimizer_status_evaluations_iterations": "exact"},
                "allocation": {"cpus_per_case": 1, "maximum_concurrent_cases": 4, "seconds_per_case": 600},
                "truth_used": False, "original_records_modified": False,
                "scope": "Production implementation compatibility; four existing diagnostic spectra and one separately reported stress spectrum; no new biological simulations or optimizer search"}
    (HERE / "manifest.json").write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")
    print("manifest_sha256", digest(HERE / "manifest.json"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--diagnostics", type=Path, required=True)
    parser.add_argument("--stress", type=Path, required=True)
    arguments = parser.parse_args()
    main(arguments.diagnostics, arguments.stress)
