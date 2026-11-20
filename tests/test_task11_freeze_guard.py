"""Synthetic provenance fixtures only; these are not validation evidence."""

import json
import sys

import pytest

from scripts import freeze_task11_multifamily_inbreeding as freeze


def evidence_fixture(tmp_path):
    source = {}
    artifacts = {}
    for paths, destination in (
        (freeze.RELEASE_SOURCES, source),
        (freeze.RELEASE_ARTIFACTS, artifacts),
    ):
        for relative in paths:
            path = tmp_path / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("unit test fixture, not a real scientific artifact\n")
            destination[relative] = freeze.sha256(path)
    records = {}
    for key in ("validation_protocol", "validation_report"):
        relative = f"test_only/{key}.json"
        path = tmp_path / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text('{"unit_test_only": true}\n')
        records[key] = {"path": relative, "sha256": freeze.sha256(path)}
    payload = {
        "schema_version": "task11-complete-validation-v1",
        "method_version": freeze.METHOD_VERSION,
        "checks": dict.fromkeys(freeze.RELEASE_CHECKS, True),
        "relationship_labels_evaluated": list(freeze.DOG_LABELS.values()),
        "confirmation_cohort_ids": ["unit-test-only-cohort"],
        "source_sha256": source,
        "artifact_sha256": artifacts,
        **records,
    }
    payload["execution_assets"] = {}
    for key in ("dog_recombination_map", "hapibd_jar", "inference_config"):
        relative = f"test_only/{key}.fixture"
        (tmp_path / relative).write_text("Synthetic unit-test execution asset")
        payload["execution_assets"][key] = {
            "path": relative,
            "sha256": freeze.sha256(tmp_path / relative),
        }
    evidence = tmp_path / "unit_test_evidence.json"
    evidence.write_text(json.dumps(payload))
    return evidence, payload


def test_freeze_without_evidence_stops_before_historical_reads(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["freeze_task11_multifamily_inbreeding.py"])

    def forbidden_read(*args, **kwargs):
        pytest.fail("Historical or final inputs were read before the evidence gate")

    monkeypatch.setattr(freeze.Path, "read_text", forbidden_read)
    monkeypatch.setattr(freeze.Path, "open", forbidden_read)
    with pytest.raises(RuntimeError, match="validation-evidence is required"):
        freeze.main()


def test_current_bound_evidence_passes_provenance_check(tmp_path):
    path, _ = evidence_fixture(tmp_path)
    result = freeze.require_release_validation(path, root=tmp_path)
    assert result["provenance_verified"] is True


@pytest.mark.parametrize(
    "change", ["source", "candidate", "model", "report", "map", "hapibd", "config"]
)
def test_modified_bound_artifact_invalidates_evidence(tmp_path, change):
    path, payload = evidence_fixture(tmp_path)
    relative = {
        "source": freeze.RELEASE_SOURCES[0],
        "candidate": freeze.RELEASE_ARTIFACTS[0],
        "model": next(
            item for item in freeze.RELEASE_ARTIFACTS if item.endswith("model.joblib")
        ),
        "report": payload["validation_report"]["path"],
        "map": payload["execution_assets"]["dog_recombination_map"]["path"],
        "hapibd": payload["execution_assets"]["hapibd_jar"]["path"],
        "config": payload["execution_assets"]["inference_config"]["path"],
    }[change]
    (tmp_path / relative).write_text("changed after validation")
    with pytest.raises(ValueError, match="hash mismatch"):
        freeze.require_release_validation(path, root=tmp_path)


@pytest.mark.parametrize(
    "change",
    [
        "method",
        "development_cohort",
        "missing_labels",
        "false_check",
        "string_check",
        "missing_binding",
    ],
)
def test_incomplete_or_development_evidence_cannot_certify_release(tmp_path, change):
    path, payload = evidence_fixture(tmp_path)
    if change == "method":
        payload["method_version"] = "task11-guarded-pair-proposal-development-v1"
    elif change == "development_cohort":
        payload["confirmation_cohort_ids"] = ["seed_2026112711"]
    elif change == "missing_labels":
        payload["relationship_labels_evaluated"] = ["siblings"]
    elif change in ("false_check", "string_check"):
        payload["checks"][freeze.RELEASE_CHECKS[0]] = (
            False if change == "false_check" else "true"
        )
    else:
        payload["artifact_sha256"].pop(freeze.RELEASE_ARTIFACTS[0])
    path.write_text(json.dumps(payload))
    with pytest.raises(ValueError):
        freeze.require_release_validation(path, root=tmp_path)
