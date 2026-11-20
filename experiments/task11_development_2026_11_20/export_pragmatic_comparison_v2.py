"""Copy complete old-model/corrected-development evidence for independent audit."""
import hashlib
import json
import shutil
import sys
import tarfile
from pathlib import Path


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


root = Path(sys.argv[1]).resolve()
result = root / "pragmatic_crossed_comparison_v2"
report = load(result / "report.json")
if report["cohorts"] != 54:
    raise ValueError("Incomplete comparison")
for path, digest in report["inputs_sha256"].items():
    if sha(path) != digest:
        raise ValueError(f"Changed input: {path}")
for name, digest in report["outputs_sha256"].items():
    if sha(result / name) != digest:
        raise ValueError(f"Changed output: {name}")
metadata = root / "pragmatic_crossed_raw_metadata_v2"
metadata.mkdir(exist_ok=False)
manifest = {}
for index in range(54):
    original = root / "results_crossed_v1" / f"development_{index:03d}"
    destination = metadata / f"development_{index:03d}"
    destination.mkdir()
    value = load(original / "result.json")
    if value["cohort"]["index"] != index or value["cohort"]["split"] != "development" or not value["completed"]:
        raise ValueError("Incorrect raw cohort identity/completion")
    for name in ("truth.json", "clusters.json", "protocol.json", "result.json"):
        digest = sha(original / name)
        if name != "result.json" and value["output_sha256"][name] != digest:
            raise ValueError("Changed raw cohort metadata")
        shutil.copyfile(original / name, destination / name)
        manifest[f"development_{index:03d}/{name}"] = digest
(metadata / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
archive = root / "pragmatic_crossed_comparison_v2_compact.tar.gz"
with tarfile.open(archive, "x:gz") as tar:
    for name in ("pragmatic_crossed_comparison_v2", "pragmatic_crossed_raw_metadata_v2",
                 "pragmatic_crossed_comparison_plan_v2.json", "crossed_development_plan_v1.json",
                 "pragmatic_context_candidate_v1.json", "evaluate_pragmatic_on_crossed_v2.py",
                 "test_pragmatic_jaccard_v2.py", "export_pragmatic_comparison_v2.py"):
        tar.add(root / name, arcname=name)
print(json.dumps({"archive": str(archive), "sha256": sha(archive), "cohorts": 54, "raw_metadata_files": len(manifest)}))
