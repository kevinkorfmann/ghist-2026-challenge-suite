"""Bind the existing context-v1 model without claiming strict qualification.

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
BASE_CANDIDATE_SHA = "354e49276bdee2643eb2b2827a397a0435c80553b7767139f7c0b07228c8a933"
CONTEXT_MODEL_SHA = "746d2715fdb875097b27b40082218a35cd5d5beb9c0dd68bef0d0aedf686b37a"
SCHEMA = "task11-pragmatic-existing-context-v1"


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


def validate_candidate(candidate, root, report, base):
    expected = {
        "schema": SCHEMA, "strict_qualification_passed": False,
        "final_candidate_selected": False, "selected_role": "context",
        "base_model": str(root / "candidate_v5/model.joblib"),
        "context_model": str(root / "context_diagnostic_v1/diagnostic_model.joblib"),
        "base_temperature": base["temperature"], "context_temperature": report["temperature"],
    }
    if any(candidate.get(key) != value for key, value in expected.items()):
        raise ValueError("Unexpected pragmatic candidate identity/model/temperature/scope")
    if (candidate.get("gates") != base["gates"] or
            candidate.get("protocol_sha256") != base["protocol_sha256"] or
            candidate["bound_inputs"].get(candidate["context_model"]) != CONTEXT_MODEL_SHA or
            candidate["bound_inputs"].get(candidate["base_model"]) != base["bound_inputs"][base["model_path"]]):
        raise ValueError("Changed strict gates, protocol or model bindings")


def initialize(root, candidate_path=None, candidate_sha=None):
    runner_path = root / "run_crossed_future_v1.py"
    verify({str(runner_path): RUNNER_SHA, str(root / "candidate_v5/candidate.json"): BASE_CANDIDATE_SHA,
            str(root / "context_diagnostic_v1/diagnostic_model.joblib"): CONTEXT_MODEL_SHA})
    spec = importlib.util.spec_from_file_location("pragmatic_verified_source_initializer", runner_path)
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    # No candidate is supplied: this verifies only the scientific source closure.
    # It does not call or weaken the strict corrected-development launch guard.
    _, protocol, _, ctx, adapter, infer = runner.initialize(root)
    base = load(root / "candidate_v5/candidate.json")
    report = load(root / "context_diagnostic_v1/report.json")
    verify(base["bound_inputs"])
    verify(report["inputs_sha256"])
    verify({str(root / "context_diagnostic_v1" / name): value for name, value in report["outputs_sha256"].items()})
    if (base["sibling_pooling"] is not False or report["parameters"] != ctx.PARAMETERS or
            report["outputs_sha256"]["diagnostic_model.joblib"] != CONTEXT_MODEL_SHA):
        raise ValueError("Unexpected existing context specification")
    candidate = None
    if candidate_path is not None:
        verify({str(candidate_path): candidate_sha})
        candidate = load(candidate_path)
        validate_candidate(candidate, root, report, base)
        verify(candidate["bound_inputs"])
        if candidate["bound_inputs"].get(str(Path(__file__).resolve())) != sha(__file__):
            raise ValueError("Unbound pragmatic initializer")
    return protocol, candidate, ctx, adapter, infer


def freeze(root):
    protocol, _, ctx, _, _ = initialize(root)
    base = load(root / "candidate_v5/candidate.json")
    report_path = root / "context_diagnostic_v1/report.json"
    report = load(report_path)
    files = [root / name for name in (
        "candidate_v5/candidate.json", "protocol_v5.json", "crossed_development_plan_v1.json",
        "context_diagnostic_v1/report.json", "context_diagnostic_source_v1/diagnose_family_context.py",
        "run_crossed_future_v1.py", "run_sparse_context_stress_v1.py", "adapter_dependency_v3.py",
        "pragmatic_context_bindings_v1.py", "prepare_pragmatic_context_v1.py", "smoke_pragmatic_context_v1.py")]
    bounds = {str(path): sha(path) for path in files}
    bounds.update(base["bound_inputs"])
    bounds.update(report["inputs_sha256"])
    bounds.update({str(root / "context_diagnostic_v1" / name): value for name, value in report["outputs_sha256"].items()})
    candidate = {"schema": SCHEMA, "frozen_utc": datetime.datetime.now(datetime.UTC).isoformat(),
        "strict_qualification_passed": False, "final_candidate_selected": False,
        "selected_role": "context", "base_model": base["model_path"],
        "context_model": str(root / "context_diagnostic_v1/diagnostic_model.joblib"),
        "base_temperature": base["temperature"], "context_temperature": report["temperature"],
        "protocol_sha256": sha(root / "protocol_v5.json"), "source_sha256": protocol["source_sha256"],
        "gates": protocol["gates"], "context_parameters": ctx.PARAMETERS, "bound_inputs": bounds,
        "historical_scores": {split: {role: value[role] for role in ("base", "context", "flips")}
                              for split, value in report["reports"].items() if split != "train"},
        "limitations": ["Development-selected fallback; no untouched confirmation passed.",
            "Exposed eight-cohort half_avuncular recall is 39/78=0.50, below the unchanged 0.60 gate.",
            "Calibration half_avuncular and avuncular_once_removed recall also fail 0.60.",
            "Old generator has ancestry/Ne confounding and sibling-rich sampling; transfer is uncertain.",
            "Model probabilities are empirical synthetic-domain classifier probabilities, not calibrated dog pedigree posteriors.",
            "One final submission is allowed; root must select the fixed candidate before final-score feedback."],
        "submission_ready": False, "submission_count": 0, "no_final_challenge_access": True}
    validate_candidate(candidate, root, report, base)
    verify(bounds)
    destination = root / "pragmatic_context_candidate_v1.json"
    with destination.open("x") as stream:
        stream.write(json.dumps(candidate, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"candidate": str(destination), "sha256": sha(destination)}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("command", choices=("freeze",))
    args = parser.parse_args()
    if not sys.flags.isolated:
        raise RuntimeError("Invoke this binding tool with python -I")
    freeze(args.root.resolve())
