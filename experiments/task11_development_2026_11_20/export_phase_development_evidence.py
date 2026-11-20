"""Reproduce frozen v4/v5 development predictions and export independent scorer input."""
import argparse
import json
from pathlib import Path
import sys

p = argparse.ArgumentParser()
p.add_argument("root", type=Path)
p.add_argument("version", choices=("v4", "v5"))
args = p.parse_args()
root, version = args.root, args.version
source = root / f"source_{version}"
sys.path[:0] = [str(source), str(source / "src"), str(source / "scripts")]
from scripts import develop_task11_pedigree as dev  # noqa: E402

protocol_path = root / f"protocol_{version}.json"
protocol = dev.load(protocol_path)
dev.require_runtime(protocol)
dev.require_hashes({str(source / path): h for path, h in protocol["source_sha256"].items()})
candidate_path = root / f"candidate_{version}/candidate.json"
candidate = dev.load(candidate_path)
dev.require_hashes(candidate["bound_inputs"])
report_path = root / f"phase_development_ablation_{version}.json"
report = dev.load(report_path)
assert report["candidate_sha256"] == dev.digest(candidate_path)
assert report["protocol_sha256"] == dev.digest(protocol_path)
model = dev.joblib.load(candidate["model_path"])
destination = root / f"raw_development_{version}"
destination.mkdir(exist_ok=False)
cohorts, bounds = [], {}
for index, expected in enumerate(report["cohorts"]):
    directory = root / "results_v3" / f"confirmation_{index:02d}"
    result, truth, rows, base, *_ = dev.dataset(directory / "result.json")
    phase_path = Path(expected["phase_features_path"])
    dev.require_hashes({str(phase_path): expected["phase_features_sha256"]})
    phase = dev.np.load(phase_path)["features"]
    predictions = dev.predict(rows, dev.np.concatenate([base, phase], axis=1), model,
                              candidate["temperature"], candidate["sibling_pooling"])
    assert dev.score(truth, predictions) == expected["phase"]
    cohorts.append({"cohort_id": f"exposed_{index:02d}", "truth": truth["pairs"],
                    "predictions": predictions, "cohort": result["cohort"]})
    bounds[str(directory / "result.json")] = dev.digest(directory / "result.json")
    bounds[str(phase_path)] = dev.digest(phase_path)
    dev.dump(destination / f"truth_{index:02d}.json", truth)
    dev.dump(destination / f"clusters_{index:02d}.json", dev.load(directory / "clusters.json"))
    dev.dump(destination / f"predictions_{index:02d}.json", predictions)
dev.dump(destination / "independent_scoring_input.json", {"labels": list(dev.LABELS), "cohorts": cohorts})
for filename, obj in (("protocol.json", protocol), ("candidate.json", candidate), ("ablation.json", report)):
    dev.dump(destination / filename, obj)
# The copied JSON is reformatted canonically by the same writer; preserve exact bindings too.
for path in (protocol_path, candidate_path, report_path, Path(__file__)):
    bounds[str(path)] = dev.digest(path)
dev.require_hashes(bounds)
dev.dump(destination / "manifest.json", {"scope": "Exposed development only; not new confirmation",
    "inputs_sha256": bounds, "outputs_sha256": {p.name: dev.digest(p) for p in sorted(destination.iterdir())},
    "source_archive_sha256": dev.digest(root / f"source_snapshot_{version}.tar.gz"),
    "protocol_sha256": dev.digest(protocol_path), "candidate_sha256": dev.digest(candidate_path),
    "model_sha256": dev.digest(Path(candidate["model_path"])), "submission_ready": False})
print(json.dumps({"version": version, "output": str(destination), "manifest_sha256": dev.digest(destination / "manifest.json")}))
