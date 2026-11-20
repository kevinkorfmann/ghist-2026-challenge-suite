"""Compact raw evidence for all216 corrected cohorts; no genotype/model copies."""
import hashlib
import json
import shutil
import sys
import tarfile
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


root = Path(sys.argv[1]).resolve()
report_path = root / "crossed_context_fit_v1/report.json"
report, plan = load(report_path), load(root / "crossed_development_plan_v1.json")
audit = load(root / "crossed_development_audit_v1/audit.json")
if audit["fit_report_sha256"] != sha(report_path):
    raise ValueError("Audit/fit report mismatch")
out = root / "crossed_all_raw_metadata_v1"
out.mkdir(exist_ok=False)
cohorts = []
for config in plan["cohorts"]:
    name = f"{config['split']}_{config['index']:03d}"
    directory = root / "results_crossed_v1" / name
    result_path = directory / "result.json"
    if report["inputs_sha256"][str(result_path)] != sha(result_path):
        raise ValueError("Bound cohort result changed")
    result = load(result_path)
    if not result["completed"] or result["cohort"] != config:
        raise ValueError("Incomplete or mismatched cohort")
    target = out / name
    target.mkdir()
    for filename in ("truth.json", "clusters.json", "protocol.json"):
        if sha(directory / filename) != result["output_sha256"][filename]:
            raise ValueError(f"Changed retained metadata: {directory / filename}")
        shutil.copyfile(directory / filename, target / filename)
    shutil.copyfile(result_path, target / "result.json")
    cohorts.append({"cohort": config, "result_sha256": sha(result_path)})
if len(cohorts) != 216 or len({c["cohort"]["seed"] for c in cohorts}) != 216:
    raise ValueError("Missing or repeated independent cohort")
manifest = {"scope": "All216 corrected training/calibration/development metadata, no outcome exclusions",
    "fit_report_sha256": sha(report_path), "plan_sha256": sha(root / "crossed_development_plan_v1.json"),
    "audit_sha256": sha(root / "crossed_development_audit_v1/audit.json"),
    "export_script_sha256": sha(__file__), "cohorts": cohorts,
    "output_sha256": {str(p.relative_to(out)): sha(p) for p in sorted(out.rglob("*.json"))}}
(out / "manifest.json").write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")
archive = root / "crossed_development_compact_v1.tar.gz"
with tarfile.open(archive, "x:gz") as tar:
    tar.add(out, arcname=out.name)
    tar.add(root / "crossed_development_audit_v1", arcname="crossed_development_audit_v1")
    for path in sorted((root / "crossed_context_fit_v1").glob("*.json")):
        if path.name in report["outputs_sha256"] and sha(path) != report["outputs_sha256"][path.name]:
            raise ValueError("Bound fit output changed")
        tar.add(path, arcname=f"crossed_context_fit_v1/{path.name}")
    for filename in ("crossed_development_plan_v1.json", "protocol_v5.json", "export_crossed_evidence_v1.py", "audit_crossed_development_v1.py"):
        tar.add(root / filename, arcname=filename)
print(json.dumps({"cohorts": len(cohorts), "manifest_sha256": sha(out / "manifest.json"),
                  "archive": str(archive), "archive_sha256": sha(archive)}))
