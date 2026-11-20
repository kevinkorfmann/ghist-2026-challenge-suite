"""Freeze one diagnostic-only retry after a cross-platform initialization assertion."""
import copy
import hashlib
import json
from pathlib import Path
import shutil


HERE = Path(__file__).resolve().parent
PARENT = HERE.parent / "production_replay"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    if (HERE / "manifest.json").exists():
        raise FileExistsError("Stress retry is already frozen")
    parent = json.loads((PARENT / "manifest.json").read_text())
    inspection = json.loads((PARENT / "stress_initialization_check.json").read_text())
    config = copy.deepcopy(parent)
    case = copy.deepcopy(parent["cases"][4])
    assert inspection["seed"] == case["seed"] == 2026112951
    case["log_starts"] = inspection["runtime_starts"]
    case["initialization"] = "Original seed; native Betty float starts frozen before retry, after inspecting only deterministic initialization"
    config.update(schema="bottleneck_production_stress_retry_v1", cases=[case], planned_cases=1,
                  planned_method_problems=2, planned_optimizer_attempts=48,
                  allocation={"cpus_per_case": 1, "maximum_concurrent_cases": 1, "seconds_per_case": 360},
                  parent_manifest_sha256=digest(PARENT / "manifest.json"),
                  initialization_inspection_sha256=digest(PARENT / "stress_initialization_check.json"),
                  scope="One diagnostic-only retry of identical production source, input, seed and optimizer budgets. No new biological or acceptance settings; not independent replication.",
                  correction="Cross-platform log/exp start rounding differs by at most8.881784197001252e-16; save fit before post-fit assertions and compare seeded runtime starts to frozen native starts.")
    shutil.copytree(PARENT / "source_snapshot", HERE / "source_snapshot")
    (HERE / "inputs").mkdir()
    shutil.copyfile(PARENT / case["input"], HERE / case["input"])
    shutil.copyfile(PARENT / "stress_initialization_check.json", HERE / "stress_initialization_check.json")
    runner = (PARENT / "run_case.py").read_text()
    runner = runner.replace('            validate_optimizer_contract(result, engine=engine,',
                            '            record["fit"] = result  # Preserve finite evidence before harness checks.\n            validate_optimizer_contract(result, engine=engine,')
    (HERE / "run_case.py").write_text(runner)
    batch = (PARENT / "run.sbatch").read_text().replace("00:10:00", "00:06:00").replace("bottleneck-production-replay-20261120", "bottleneck-production-stress-retry-20261120")
    (HERE / "run.sbatch").write_text(batch)
    source_hashes = {name: expected for name, expected in parent["source_sha256"].items() if name.startswith("source_snapshot/")}
    source_hashes.update({name: digest(HERE / name) for name in ("prepare.py", "run_case.py", "run.sbatch", "stress_initialization_check.json")})
    config["source_sha256"] = source_hashes
    (HERE / "manifest.json").write_text(json.dumps(config, indent=2, allow_nan=False) + "\n")
    print(digest(HERE / "manifest.json"))


if __name__ == "__main__":
    main()
