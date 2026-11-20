"""Bind the single expanded model to prospective stress only after fixed gates pass."""
import datetime
import hashlib
import json
import sys
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


root = Path(sys.argv[1]).resolve()
for name in ("results_sparse_context_v1", "predictions_sparse_context_v1", "sparse_context_candidate_v1.json"):
    if (root / name).exists():
        raise FileExistsError(f"Candidate must be frozen before stress: {name}")
fixed = {"sparse_context_stress_design_v1.json": "2b4f32a1f8b1083e119ebde8351c16921028f5ce79f4ec4b13c5b4c969897b0c",
    "sparse_context_stress_launch_guard_v1.json": "99e66511db306365c7611b35820d43c5f1130669093d8f589330d42a432e4587",
    "sparse_context_generation_protocol_v1.json": "81310ce40c099fb66243e0043ed4479e84ebc1e7faf6aff39c2d4f7410d92483",
    "protocol_v5.json": "ca55f4dd692e0c1aa41704dcb3901b30b129437b1e0905211ceb84f00bfd0d1c",
    "source_snapshot_v5.tar.gz": "77375babb7fdac32b16ed7e487267f1fb38b29e6f7d291bcafaa4bfc8a86f558",
    "expanded_training_plan_v2.json": "0d1c89dc3c6f8925c0c5d88653cf555e28c44343d2f0c98509f9e08fb5f98da8",
    "context_diagnostic_source_v1/diagnose_family_context.py": "4b347c189e0a5799c99e158715c236934485e618db44cb231e84a27535fd396b"}
for name, expected in fixed.items():
    if sha(root / name) != expected:
        raise ValueError(f"Prospective binding changed: {name}")
audit_path = root / "expanded_development_audit_v1/audit.json"
audit = load(audit_path)
if not audit["sparse_stress_launch_allowed"]:
    raise ValueError("Expanded development gates did not pass")
if not audit["summaries"]["exposed_context"]["all_nine_gates_passed"]:
    raise ValueError("Exposed development gates did not pass")
calibration = audit["summaries"]["calibration_context"]["checks"]
if not calibration["each_label_recall"] or not calibration["calibration"]:
    raise ValueError("Calibration launch gates did not pass")
report_path = root / "expanded_context_fit_v2/report.json"
report = load(report_path)
if sha(report_path) != audit["fit_report_sha256"]:
    raise ValueError("Audited fit report changed")
paths = [root / name for name in fixed]
paths.extend([audit_path, report_path, Path(__file__).resolve(), root / "run_sparse_context_stress_v1.py",
              root / "audit_expanded_development_v1.py", root / "sparse_context_stress_v1.sbatch",
              root / "sparse_context_summary_v1.sbatch"])
paths.extend(root / "expanded_context_fit_v2" / name for name in ("base_model.joblib", "context_model.joblib"))
for role in ("base", "context"):
    name = f"{role}_model.joblib"
    if sha(root / "expanded_context_fit_v2" / name) != report["outputs_sha256"][name]:
        raise ValueError("Expanded model changed")
bound = {str(path): sha(path) for path in paths}
for directory, hashes in ((root / "expanded_development_audit_v1", audit["output_sha256"]),
                          (root / "expanded_context_fit_v2", report["outputs_sha256"])):
    for filename, expected in hashes.items():
        path = directory / filename
        if sha(path) != expected:
            raise ValueError(f"Retained evidence changed: {path}")
        bound[str(path)] = expected
candidate = {"schema": "task11-frozen-expanded-context-stress-v1",
    "frozen_before_stress": True, "frozen_at_utc": datetime.datetime.now(datetime.UTC).isoformat(),
    "runner_sha256": sha(root / "run_sparse_context_stress_v1.py"),
    "base_model": str(root / "expanded_context_fit_v2/base_model.joblib"),
    "context_model": str(root / "expanded_context_fit_v2/context_model.joblib"),
    "base_temperature": report["base_temperature"], "context_temperature": report["context_temperature"],
    "bound_inputs": bound, "preservation_audit": audit["preservation_audit"],
    "gates": audit["gates"], "scope": "Single fixed expanded context model; additional sparse stress before reserved confirmation",
    "submission_ready": False, "final_challenge_access": False, "refitting_allowed": False}
path = root / "sparse_context_candidate_v1.json"
with path.open("x") as destination:
    json.dump(candidate, destination, indent=2)
    destination.write("\n")
print(json.dumps({"candidate_sha256": sha(path), "base_temperature": candidate["base_temperature"],
                  "context_temperature": candidate["context_temperature"], "bound_inputs": len(bound)}))
