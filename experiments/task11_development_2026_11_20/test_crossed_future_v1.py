"""Meaningful fail-closed gates for future execution; no simulated cohorts."""
import copy
import importlib.util
from pathlib import Path

import pytest

DIRECTORY = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("future_runner_test", DIRECTORY / "run_crossed_future_v1.py")
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


def development():
    return {"input_design_sha256": runner.PLAN_SHA, "fresh_stress_launch_allowed": True,
        "scientific_development_gates_passed": True, "within_decision_budget": True,
        "summaries": {"development_context": {"checks": dict.fromkeys(runner.CHECKS, True), "all_nine_gates_passed": True},
            "calibration_context": {"checks": {"each_label_recall": True, "calibration": True, "bootstrap_lower": False}}}}


def test_calibration_bootstrap_is_not_an_added_gate():
    runner.require_development(development())


@pytest.mark.parametrize("gate", sorted(runner.CHECKS))
def test_any_development_gate_failure_prevents_future(gate):
    value = development()
    value["summaries"]["development_context"]["checks"][gate] = False
    with pytest.raises(ValueError, match="launch gates failed"):
        runner.require_development(value)


@pytest.mark.parametrize("gate", ["each_label_recall", "calibration"])
def test_calibration_launch_gates(gate):
    value = development()
    value["summaries"]["calibration_context"]["checks"][gate] = False
    with pytest.raises(ValueError):
        runner.require_development(value)


def test_unknown_or_missing_gate_rejected():
    value = development()
    value["summaries"]["development_context"]["checks"]["extra"] = True
    with pytest.raises(ValueError):
        runner.require_development(value)
    del value["summaries"]["development_context"]["checks"]["extra"]
    del value["summaries"]["development_context"]["checks"]["accuracy"]
    with pytest.raises(ValueError):
        runner.require_development(value)


@pytest.mark.parametrize("field", ["fresh_stress_launch_allowed", "scientific_development_gates_passed", "within_decision_budget"])
def test_failed_or_expired_audit_rejected(field):
    value = development()
    value[field] = False
    with pytest.raises(ValueError):
        runner.require_development(value)


def test_stress_is_bound_to_candidate_and_all_eight_cohorts():
    value = {"split": "stress", "candidate_sha256": "a" * 64, "plan_sha256": runner.PLAN_SHA,
        "all_original_gates_passed": True, "checks": dict.fromkeys(runner.CHECKS, True), "cohorts": list(range(8))}
    runner.require_stress(value, "a" * 64)
    for patch in ({"candidate_sha256": "b" * 64}, {"split": "confirmation"},
                  {"cohorts": list(range(7))}, {"plan_sha256": "b" * 64}):
        with pytest.raises(ValueError):
            runner.require_stress({**value, **patch}, "a" * 64)
    failed = copy.deepcopy(value)
    failed["checks"]["each_label_recall"] = False
    with pytest.raises(ValueError):
        runner.require_stress(failed, "a" * 64)


def test_launch_requires_slurm_before_any_generation(monkeypatch, tmp_path):
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    with pytest.raises(RuntimeError, match="Slurm"):
        runner.check_launch(tmp_path, "stress", {}, "a" * 64, None)


def test_launch_deadline_cannot_be_bypassed(monkeypatch, tmp_path):
    monkeypatch.setenv("SLURM_JOB_ID", "test")
    plan = {"budget": {"latest_stress_confirmation_start_utc": "2000-01-01T00:00:00+00:00"}}
    with pytest.raises(RuntimeError, match="deadline"):
        runner.check_launch(tmp_path, "stress", plan, "a" * 64, None)


def test_hash_failure_prevents_import(tmp_path):
    path = tmp_path / "input.py"
    path.write_text("raise RuntimeError('must not execute')")
    with pytest.raises(ValueError, match="changed"):
        runner.verify({str(path): "0" * 64})


def test_full_denominator_gate_rejects_missing_and_extra_pairs():
    plan = runner.load(DIRECTORY / "crossed_development_plan_v1.json")
    perfect = {"accuracy": 1., "balanced_accuracy": 1., "recall": {label: 1. for label in plan["labels"]},
        "cohort_bootstrap_balanced_interval95": [.9, 1.], "missing_true_pairs": 1,
        "extra_cross_family_pairs": 1, "ece_10_equal_width": 0.}
    checks = runner.gate_checks(perfect, 1., 0, plan["gates"])
    assert checks["missing_true_pairs"] is False
    assert checks["extra_pairs"] is False
    assert checks["accuracy"] is True
