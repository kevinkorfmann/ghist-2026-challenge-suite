"""Retain all64 training pedigrees/partitions without copying genotype payloads."""
import hashlib
import json
import shutil
import sys
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


root = Path(sys.argv[1]).resolve()
report_path = root / "expanded_context_fit_v2/report.json"
report = json.loads(report_path.read_text())
out = root / "expanded_training_raw_metadata_v1"
out.mkdir(exist_ok=False)
cohorts = []
for index in range(64):
    directory = root / ("results_v5" if index < 16 else "results_expanded_v1") / f"train_{index:02d}"
    result_path = directory / "result.json"
    if report["inputs_sha256"][str(result_path)] != sha(result_path):
        raise ValueError("Training result changed")
    result = json.loads(result_path.read_text())
    if not result["completed"]:
        raise ValueError("Incomplete training cohort")
    target = out / f"train_{index:02d}"
    target.mkdir()
    for name in ("truth.json", "clusters.json", "protocol.json"):
        if sha(directory / name) != result["output_sha256"][name]:
            raise ValueError(f"Training raw metadata changed: {directory / name}")
        shutil.copyfile(directory / name, target / name)
    shutil.copyfile(result_path, target / "result.json")
    cohorts.append({"cohort": result["cohort"], "result_sha256": sha(result_path)})
if len({cohort["cohort"]["seed"] for cohort in cohorts}) != 64:
    raise ValueError("Repeated independent training seed")
manifest = {"scope": "All64 independent training raw metadata, no outcome exclusions",
    "fit_report_sha256": sha(report_path), "export_script_sha256": sha(__file__), "cohorts": cohorts,
    "output_sha256": {str(path.relative_to(out)): sha(path) for path in sorted(out.rglob("*.json"))}}
with (out / "manifest.json").open("x") as destination:
    json.dump(manifest, destination, indent=2)
    destination.write("\n")
print(json.dumps({"cohorts": len(cohorts), "files": len(manifest["output_sha256"]), "manifest_sha256": sha(out / "manifest.json")}))
