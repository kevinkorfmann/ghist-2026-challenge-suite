"""Read small generation metadata only; no genotypes, model fitting or scoring."""
import hashlib
import json
import sys
from pathlib import Path

root = Path(sys.argv[1])
expected = hashlib.sha256((root / "expanded_training_protocol_v1.json").read_bytes()).hexdigest()
results = []
for path in sorted((root / "results_expanded_v1").glob("train_*/result.json")):
    result = json.loads(path.read_text())
    protocol = json.loads(path.with_name("protocol.json").read_text())
    assert protocol["master_sha256"] == expected
    results.append(result)
successful = [r for r in results if r["completed"]]
print(json.dumps({"completed": len(successful), "failed": [r.get("failure") for r in results if not r["completed"]],
    "simple_completed": sum(not r["cohort"]["complex_topology"] for r in successful),
    "mixed_completed": sum(r["cohort"]["complex_topology"] for r in successful),
    "exact_families": sum(r["partition"]["exactly_recovered_families"] for r in successful),
    "total_families": len(successful)*8, "scope": "Generation metadata only; full hashes rechecked by fit"}))
