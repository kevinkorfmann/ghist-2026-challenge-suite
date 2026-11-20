"""Bind the corrected fixed base/context model without claiming strict qualification.

This separate schema cannot pass the frozen prospective qualification runner.
It preserves the failed utility gates and never chooses a final submission.
"""
import argparse
import datetime
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

RUNNER_SHA = "49ab0f17b3948fa87e416682f16cf2bc4a785424f7a1328ea95c3dfefbcedfa2"
SCHEMA = "task11-pragmatic-corrected-fixed-model-v1"


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def verify(bound):
    if not isinstance(bound, dict) or not bound:
        raise ValueError("Missing hash bindings")
    for path, expected in bound.items():
        if sha(path) != expected:
            raise ValueError(f"Bound input changed: {path}")


def validate_candidate(candidate, root, report, protocol):
    if candidate.get("selected_role") not in ("base", "context"):
        raise ValueError("Unknown fixed model role")
    expected = {"schema": SCHEMA, "strict_qualification_passed": False,
        "final_candidate_selected": False, "protocol_sha256": sha(root / "protocol_v5.json"),
        "gates": protocol["gates"]}
    for role in ("base", "context"):
        model = root / "crossed_context_fit_v1" / f"{role}_model.joblib"
        expected[f"{role}_model"] = str(model)
        expected[f"{role}_temperature"] = report[f"{role}_temperature"]
        if candidate["bound_inputs"].get(str(model)) != report["outputs_sha256"][model.name]:
            raise ValueError("Changed fitted model binding")
    if any(candidate.get(key) != value for key, value in expected.items()):
        raise ValueError("Changed candidate identity, model, temperature or scientific scope")


def initialize(root, candidate_path=None, candidate_sha=None):
    runner_path = root / "run_crossed_future_v1.py"
    verify({str(runner_path): RUNNER_SHA})
    spec = importlib.util.spec_from_file_location("pragmatic_verified_source_initializer", runner_path)
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    # No candidate is supplied: this verifies only the scientific source closure.
    # It does not call or weaken the strict corrected-development launch guard.
    _, protocol, _, ctx, adapter, infer = runner.initialize(root)
    report = load(root / "crossed_context_fit_v1/report.json")
    verify(report["inputs_sha256"])
    verify({str(root / "crossed_context_fit_v1" / name): value for name, value in report["outputs_sha256"].items()})
    audit = load(root / "crossed_development_audit_v1/audit.json")
    if audit["fit_report_sha256"] != sha(root / "crossed_context_fit_v1/report.json"):
        raise ValueError("Raw development audit differs from fit")
    if report["training_audit"]["cohorts"] != 108 or report["training_audit"]["families"] != 864:
        raise ValueError("Incomplete corrected training")
    for split in ("calibration", "development"):
        for role in ("base", "context"):
            if len(load(root / "crossed_context_fit_v1" / f"scoring_{split}_{role}.json")["cohorts"]) != 54:
                raise ValueError("Incomplete corrected calibration/development")
    candidate = None
    if candidate_path is not None:
        verify({str(candidate_path): candidate_sha})
        candidate = load(candidate_path)
        validate_candidate(candidate, root, report, protocol)
        verify(candidate["bound_inputs"])
        if candidate["bound_inputs"].get(str(Path(__file__).resolve())) != sha(__file__):
            raise ValueError("Unbound pragmatic initializer")
    return protocol, candidate, ctx, adapter, infer


def freeze(root, role):
    protocol, _, ctx, _, _ = initialize(root)
    report_path = root / "crossed_context_fit_v1/report.json"
    report = load(report_path)
    audit_path = root / "crossed_development_audit_v1/audit.json"
    audit = load(audit_path)
    files = [root / name for name in (
        "protocol_v5.json", "crossed_development_plan_v1.json", "crossed_context_fit_v1/report.json",
        "crossed_development_audit_v1/audit.json", "context_diagnostic_source_v1/diagnose_family_context.py",
        "run_crossed_future_v1.py", "run_sparse_context_stress_v1.py", "adapter_dependency_v3.py",
        "pragmatic_corrected_bindings_v1.py", "prepare_pragmatic_corrected_v1.py", "smoke_pragmatic_corrected_v1.py")]
    bounds = {str(path): sha(path) for path in files}
    bounds.update(report["inputs_sha256"])
    bounds.update({str(root / "crossed_context_fit_v1" / name): value for name, value in report["outputs_sha256"].items()})
    candidate = {"schema": SCHEMA, "frozen_utc": datetime.datetime.now(datetime.UTC).isoformat(),
        "strict_qualification_passed": False, "final_candidate_selected": False,
        "selected_role": role, "protocol_sha256": sha(root / "protocol_v5.json"),
        "source_sha256": protocol["source_sha256"], "gates": protocol["gates"],
        "context_parameters": ctx.PARAMETERS, "bound_inputs": bounds,
        "strict_development_audit": audit,
        "limitations": ["Packaging a fixed model does not select the final submission or certify independent confirmation.",
            "Any failed strict development gates are retained verbatim in the bound raw audit.",
            "Corrected finite synthetic templates cover1–4 member branches and varied inbreeding/populations, not all dog pedigrees.",
            "True phase, MCAR missingness and coalescent founder assumptions limit transfer.",
            "Classifier probabilities are synthetic-domain estimates, not validated dog pedigree posteriors.",
            "Final single-submission selection remains root-owned and must not use final-score feedback."],
        "submission_ready": False, "submission_count": 0, "no_final_challenge_access": True}
    for model_role in ("base", "context"):
        candidate[f"{model_role}_model"] = str(root / "crossed_context_fit_v1" / f"{model_role}_model.joblib")
        candidate[f"{model_role}_temperature"] = report[f"{model_role}_temperature"]
    validate_candidate(candidate, root, report, protocol)
    verify(bounds)
    destination = root / f"pragmatic_corrected_{role}_candidate_v1.json"
    with destination.open("x") as stream:
        stream.write(json.dumps(candidate, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"candidate": str(destination), "sha256": sha(destination)}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("command", choices=("freeze",))
    parser.add_argument("--role", choices=("base", "context"), required=True)
    args = parser.parse_args()
    if not sys.flags.isolated:
        raise RuntimeError("Invoke this binding tool with python -I")
    freeze(args.root.resolve(), args.role)
