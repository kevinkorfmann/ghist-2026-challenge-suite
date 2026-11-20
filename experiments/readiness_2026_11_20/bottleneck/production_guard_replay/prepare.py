"""Freeze the duplicate-start guard patch on existing inputs; no truth or new biology."""
import copy
import hashlib
import json
from pathlib import Path
import shutil

from ghist_bottleneck.optimization import optimizer_contract


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
PARENT = HERE.parent / "production_replay"
STRESS = HERE.parent / "production_stress_retry"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    if (HERE / "manifest.json").exists():
        raise FileExistsError("Final guard replay already frozen")
    parent = json.loads((PARENT / "manifest.json").read_text())
    stress_parent = json.loads((STRESS / "manifest.json").read_text())
    config = copy.deepcopy(parent)
    stress = copy.deepcopy(stress_parent["cases"][0])
    stress["scope"] = "separate_stress"
    stress["expected_attempts"] = {}
    for engine in stress["methods"]:
        path = STRESS / "results" / stress["name"] / f"{engine}.json"
        result = json.loads(path.read_text())
        if result["state"] != "COMPLETED":
            raise ValueError("Prior stress retry must be fully retained")
        stress["reference_sha256"][f"production_stress_retry/results/{stress['name']}/{engine}.json"] = digest(path)
        stress["expected_attempts"][engine] = [{"arm": row["arm"], "restart": row["restart"], "parameters": row["parameters"],
            "objective": row["objective"], "optimizer": {"status": row["optimizer"]["warnflag"],
            "evaluations": row["optimizer"]["evaluations"], "iterations": row["optimizer"]["iterations"]}}
            for row in result["fit"]["optimizer_attempts"]]
    config["cases"] = copy.deepcopy(parent["cases"][:4]) + [stress]
    config.update(schema="bottleneck_production_guard_replay_v1", parent_manifest_sha256=digest(PARENT / "manifest.json"),
                  stress_retry_manifest_sha256=digest(STRESS / "manifest.json"),
                  scope="Source-current duplicate-start guard verification on identical prior inputs/starts. Four reference spectra and one separately reported stress spectrum; not independent biological confirmation.")
    (HERE / "inputs").mkdir()
    for case in config["cases"]:
        shutil.copyfile(PARENT / case["input"], HERE / case["input"])
    source_hashes = {}
    for path in sorted((PARENT / "source_snapshot/src").rglob("*.py")):
        relative = path.relative_to(PARENT / "source_snapshot/src")
        target = HERE / "source_snapshot/src" / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / "src" / relative, target)
        source_hashes[str(target.relative_to(HERE))] = digest(target)
    contracts = {engine: optimizer_contract(engine) for engine in ("moments", "dadi")}
    for contract in contracts.values():
        contract["versions"] = {name: config["runtime_versions"][name] for name in contract["versions"]}
    config["optimizer_contracts"] = contracts
    shutil.copyfile(STRESS / "run_case.py", HERE / "run_case.py")
    batch = (PARENT / "run.sbatch").read_text().replace("bottleneck-production-replay-20261120", "bottleneck-production-guard-replay-20261120")
    (HERE / "run.sbatch").write_text(batch)
    config["source_sha256"] = {**source_hashes, **{name: digest(HERE / name) for name in ("prepare.py", "run_case.py", "run.sbatch")}}
    (HERE / "manifest.json").write_text(json.dumps(config, indent=2, allow_nan=False) + "\n")
    print(digest(HERE / "manifest.json"))


if __name__ == "__main__":
    main()
