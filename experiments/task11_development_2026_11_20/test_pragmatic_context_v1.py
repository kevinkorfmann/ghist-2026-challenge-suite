"""Fail-closed identity and scope checks for the separate pragmatic package."""
import importlib.util
from pathlib import Path

import pytest

spec = importlib.util.spec_from_file_location("pragmatic_bindings_test", Path(__file__).with_name("pragmatic_context_bindings_v1.py"))
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def fixture():
    root = Path("/synthetic")
    base = {"temperature": .75, "gates": {"each_label_recall_min": .6}, "protocol_sha256": "protocol",
            "model_path": str(root / "candidate_v5/model.joblib"), "bound_inputs": {str(root / "candidate_v5/model.joblib"): "base"}}
    report = {"temperature": .75}
    candidate = {"schema": mod.SCHEMA, "strict_qualification_passed": False, "final_candidate_selected": False,
        "selected_role": "context", "base_model": base["model_path"],
        "context_model": str(root / "context_diagnostic_v1/diagnostic_model.joblib"),
        "base_temperature": .75, "context_temperature": .75, "gates": base["gates"],
        "protocol_sha256": "protocol", "bound_inputs": {base["model_path"]: "base",
            str(root / "context_diagnostic_v1/diagnostic_model.joblib"): mod.CONTEXT_MODEL_SHA}}
    return root, base, report, candidate


def test_valid_existing_unqualified_identity():
    root, base, report, candidate = fixture()
    mod.validate_candidate(candidate, root, report, base)


@pytest.mark.parametrize("field,value", [
    ("schema", "task11-fixed-crossed-context-candidate-v1"),
    ("strict_qualification_passed", True), ("final_candidate_selected", True),
    ("selected_role", "base"), ("base_model", "/other"), ("context_model", "/other"),
    ("base_temperature", 1), ("context_temperature", 1),
    ("protocol_sha256", "other"), ("gates", {"each_label_recall_min": .5})])
def test_reject_changed_identity_or_scientific_claim(field, value):
    root, base, report, candidate = fixture()
    candidate[field] = value
    with pytest.raises(ValueError):
        mod.validate_candidate(candidate, root, report, base)


@pytest.mark.parametrize("role", ["base", "context"])
def test_reject_changed_model_hash(role):
    root, base, report, candidate = fixture()
    candidate["bound_inputs"][candidate[f"{role}_model"]] = "different"
    with pytest.raises(ValueError):
        mod.validate_candidate(candidate, root, report, base)


def test_verify_rejects_mutated_bound_bytes(tmp_path):
    path = tmp_path / "source.py"
    path.write_text("pass\n")
    bound = {str(path): mod.sha(path)}
    mod.verify(bound)
    path.write_text("raise RuntimeError()\n")
    with pytest.raises(ValueError):
        mod.verify(bound)


def test_missing_bindings_fail():
    with pytest.raises(ValueError):
        mod.verify({})
