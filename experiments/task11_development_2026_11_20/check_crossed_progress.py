"""Read-only metadata progress for the frozen crossed generation array."""
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

root = Path(sys.argv[1])
plan_path = root / "crossed_development_plan_v1.json"
plan = json.loads(plan_path.read_text())
plan_hash = hashlib.sha256(plan_path.read_bytes()).hexdigest()
completed, failed, exact, kinds_checked = Counter(), [], 0, 0
for config in plan["cohorts"]:
    directory = root / "results_crossed_v1" / f"{config['split']}_{config['index']:03d}"
    path = directory / "result.json"
    if not path.exists():
        continue
    result = json.loads(path.read_text())
    assert result["cohort"] == config
    assert json.loads((directory / "protocol.json").read_text())["master_sha256"] == plan_hash
    if not result["completed"]:
        failed.append({"cohort": config["global_index"], "failure": result.get("failure")})
        continue
    truth = json.loads((directory / "truth.json").read_text())
    assert [p["kind"] for p in truth["pedigrees"]] == config["ancestry_kinds"]
    kinds_checked += 8
    completed[config["split"]] += 1
    exact += result["partition"]["exactly_recovered_families"]
print(json.dumps({"completed": sum(completed.values()), "by_split": dict(completed), "failed": failed,
                  "exact_families": exact, "families_checked_for_ancestry": kinds_checked,
                  "scope": "Read-only generation metadata; full retained hashes checked by fit"}))
