"""Fail-closed identity and scope checks for the separate pragmatic package."""
import importlib.util
from pathlib import Path

import pytest

spec = importlib.util.spec_from_file_location("pragmatic_bindings_test", Path(__file__).with_name("pragmatic_corrected_bindings_v1.py"))
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def fixture(tmp_path):
    root = tmp_path
    (root / "protocol_v5.json").write_text("{}")
    protocol = {"gates": {"each_label_recall_min": .6}}
    report = {"base_temperature": .75, "context_temperature": 1.,
        "outputs_sha256": {"base_model.joblib": "base", "context_model.joblib": "context"}}
    candidate = {"schema": mod.SCHEMA, "strict_qualification_passed": False, "final_candidate_selected": False,
        "selected_role": "context", "protocol_sha256": mod.sha(root / "protocol_v5.json"),
        "gates": protocol["gates"], "bound_inputs": {}}
    for role in ("base", "context"):
        candidate[f"{role}_model"] = str(root / "crossed_context_fit_v1" / f"{role}_model.joblib")
        candidate[f"{role}_temperature"] = report[f"{role}_temperature"]
        candidate["bound_inputs"][candidate[f"{role}_model"]] = role
    return root, protocol, report, candidate


def test_valid_existing_unqualified_identity(tmp_path):
    root, base, report, candidate = fixture(tmp_path)
    mod.validate_candidate(candidate, root, report, base)


@pytest.mark.parametrize("field,value", [
    ("schema", "task11-fixed-crossed-context-candidate-v1"),
    ("strict_qualification_passed", True), ("final_candidate_selected", True),
    ("selected_role", "unknown"), ("base_model", "/other"), ("context_model", "/other"),
    ("base_temperature", 1), ("context_temperature", .75),
    ("protocol_sha256", "other"), ("gates", {"each_label_recall_min": .5})])
def test_reject_changed_identity_or_scientific_claim(field, value, tmp_path):
    root, base, report, candidate = fixture(tmp_path)
    candidate[field] = value
    with pytest.raises(ValueError):
        mod.validate_candidate(candidate, root, report, base)


@pytest.mark.parametrize("role", ["base", "context"])
def test_reject_changed_model_hash(role, tmp_path):
    root, base, report, candidate = fixture(tmp_path)
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


def test_base_role_is_an_explicit_fixed_alternative(tmp_path):
    root, protocol, report, candidate = fixture(tmp_path)
    candidate["selected_role"] = "base"
    mod.validate_candidate(candidate, root, report, protocol)
