from __future__ import annotations

import csv
import importlib.metadata
import json
from pathlib import Path
import platform
import subprocess
import sys

import joblib
import numpy as np
import pytest
from sklearn.ensemble import RandomForestClassifier

from ghist_relatedness.ibd import matrix, pair_features, read_sex_map
from ghist_relatedness.labels import IBD_LABEL_MEIOSES
from ghist_relatedness.task9_release import input_pairs
from scripts.prepare_task9_submission_package import (
    EVIDENCE, RELEASE, REQUIREMENTS, SOURCES, digest, prepare_package, preflight, release_contract,
    run_final, verify_package, write_json,
)


ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def package(tmp_path):
    destination = tmp_path / "package"
    prepare_package(ROOT, ROOT / RELEASE, destination)
    return destination


def test_offline_package_leaves_final_fields_pending_without_model(package):
    manifest = verify_package(package)
    assert manifest["final_input"] == {
        "entity": "syn75202631", "version": None, "sha256": None, "path": None,
    }
    assert manifest["final_prediction"]["path"] is None
    assert manifest["confirmation_pairs"] == 20023
    assert not list(package.rglob("*.joblib"))
    assert not list(package.rglob("*.csv"))


@pytest.mark.parametrize("change", [
    {"schema": "historical.freeze"}, {"status": "CANDIDATE_NOT_FROZEN"},
    {"model": {"path": "artifacts/task9_ibd_segments/model.joblib"}},
    {"gates": {"accuracy": False}}, {"release_id": "../old"},
    {"gates": {"accuracy": True}}, {"feature_contract_sha256": {}},
])
def test_historical_or_default_release_is_refused(change):
    release = json.loads((ROOT / RELEASE).read_text())
    release.update(change)
    with pytest.raises(ValueError, match="historical/default"):
        release_contract(release)


def test_pinned_predictor_tampering_is_rejected(package):
    (package / "scripts/predict_task9_current.py").write_text("pass\n")
    with pytest.raises(ValueError, match="checksum mismatch"):
        verify_package(package)


@pytest.mark.parametrize("field,value", [
    ("runtime", {"python": "0.0.0", "packages": {}}),
    ("external_artifacts", {}), ("training_pairs", 1), ("confirmation_pairs", 1),
    ("submission_count", 1), ("operational_passed", False),
    ("final_input", {"entity": "syn75202631", "version": 1}),
    ("status", "LOCAL_FINAL_FILE_VALIDATED_NOT_SUBMITTED"),
])
def test_manifest_metadata_must_match_pinned_files(package, field, value):
    manifest = json.loads((package / "package.json").read_text())
    manifest[field] = value
    write_json(package / "package.json", manifest)
    with pytest.raises(ValueError, match="metadata differs"):
        verify_package(package)


def test_external_manifest_pin_detects_manifest_replacement(package):
    with pytest.raises(ValueError, match="externally pinned checksum"):
        verify_package(package, "0" * 64)


def test_rehashed_evidence_with_different_model_is_rejected(package):
    evidence = package / EVIDENCE[4]
    report = json.loads(evidence.read_text())
    report["model_sha256"] = "0" * 64
    write_json(evidence, report)
    manifest = json.loads((package / "package.json").read_text())
    manifest["files_sha256"][EVIDENCE[4]] = digest(evidence)
    write_json(package / "package.json", manifest)
    with pytest.raises(ValueError, match="different model"):
        verify_package(package)


def test_rehashed_nonfinite_operational_metric_is_rejected(package):
    evidence = package / EVIDENCE[4]
    report = json.loads(evidence.read_text())
    report["datasets"][0]["cases"]["official_csv"]["accuracy"] = float("nan")
    write_json(evidence, report)
    manifest = json.loads((package / "package.json").read_text())
    manifest["files_sha256"][EVIDENCE[4]] = digest(evidence)
    write_json(package / "package.json", manifest)
    with pytest.raises(ValueError, match="identity/coverage/accuracy mismatch"):
        verify_package(package)


def test_package_symlink_is_refused_even_when_content_matches(package):
    source = package / SOURCES[1]
    source.unlink()
    source.symlink_to(ROOT / SOURCES[1])
    with pytest.raises(ValueError, match="symlink refused"):
        verify_package(package)


@pytest.mark.parametrize("name", ["src/unpinned.py", "src/ghist_relatedness/__pycache__/ibd.cpython-312.pyc"])
def test_unpinned_import_files_are_refused(package, name):
    path = package / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"untrusted cache or import")
    with pytest.raises(ValueError, match="unpinned bytecode/imports"):
        verify_package(package)


def test_manifest_symlink_is_refused(package, tmp_path):
    manifest = package / "package.json"
    outside = tmp_path / "manifest.json"
    manifest.rename(outside)
    manifest.symlink_to(outside)
    with pytest.raises(ValueError, match="symlink refused"):
        verify_package(package)


@pytest.mark.parametrize("action", ["prepare", "run"])
def test_dangling_output_symlink_does_not_create_target(package, tmp_path, action):
    destination = tmp_path / "dangling"
    target = tmp_path / "must-not-be-created"
    destination.symlink_to(target, target_is_directory=True)
    with pytest.raises(ValueError, match="already exists"):
        if action == "prepare":
            prepare_package(ROOT, ROOT / RELEASE, destination)
        else:
            run_final(package, "missing", "missing", "missing", "missing", destination)
    assert not target.exists()


def test_existing_package_is_not_replaced(package):
    with pytest.raises(ValueError, match="already exists"):
        prepare_package(ROOT, ROOT / RELEASE, package)


def test_mismatched_evidence_is_rejected(tmp_path):
    release = json.loads((ROOT / RELEASE).read_text())
    release["evidence_sha256"] = {"not the evidence": "0" * 64}
    path = tmp_path / "release.json"
    write_json(path, release)
    with pytest.raises(ValueError, match="evidence differs"):
        prepare_package(ROOT, path, tmp_path / "package")


def test_stale_run_directory_is_refused_before_reading_input(package, tmp_path):
    stale = tmp_path / "stale"
    stale.mkdir()
    (stale / "task9_final_candidate.csv").write_text("historical candidate")
    with pytest.raises(ValueError, match="historical candidates"):
        run_final(package, "missing", "missing", "missing", "missing", stale)


@pytest.mark.parametrize("changes", [
    {"entity": "syn75202632"}, {"version": True}, {"version": 0}, {"sha256": "0" * 64},
    {"entity_id": "syn75202632"},
])
def test_input_provenance_rejected_before_model_access(package, tmp_path, changes):
    observed = tmp_path / "synthetic.csv"
    observed.write_text("synthetic fixture, not challenge data\n")
    provenance = {"entity": "syn75202631", "version": 1, "sha256": digest(observed)}
    provenance.update(changes)
    path = tmp_path / "provenance.json"
    write_json(path, provenance)
    with pytest.raises(ValueError, match="entity/version/checksum"):
        run_final(package, "missing", "missing", observed, path, tmp_path / "run")
    assert not (tmp_path / "run").exists()


def test_incompatible_runtime_is_rejected_before_model_access(package, monkeypatch):
    monkeypatch.setattr(platform, "python_version", lambda: "0.0.0")
    with pytest.raises(ValueError, match="recorded release runtime"):
        preflight(package, "missing", "missing")


def test_documented_feature_guard_limits_remain_explicit(tmp_path):
    genetic_map = tmp_path / "synthetic.map"
    genetic_map.write_text("1 0 0 0\n1 20 5 15\n")
    row = {"indID1": "synthetic_A", "indID2": "synthetic_B", "chr": "1", "start_bp": "0",
           "end_bp": "20", "type": "IBD1", "start_cM": "0", "end_cM": "10", "length": "10"}
    clamped = {**row, "end_bp": "200"}
    reference = pair_features([row], 3350.933402, read_sex_map(genetic_map))
    assert input_pairs([clamped]) == {("synthetic_A", "synthetic_B")}
    assert pair_features([clamped], 3350.933402, read_sex_map(genetic_map)) == reference
    fractional_zero = {**row, "start_bp": "0.5", "end_bp": "0.5", "end_cM": "0", "length": "0"}
    overlap = {**row, "start_bp": "1", "start_cM": "1", "length": "9"}
    assert input_pairs([row, fractional_zero, overlap]) == {("synthetic_A", "synthetic_B")}


@pytest.fixture
def synthetic_run(package, tmp_path):
    # All inputs/model below are generated here. No current challenge data is used.
    row = {"indID1": "synthetic_A", "indID2": "synthetic_B", "chr": "1", "start_bp": "0",
           "end_bp": "20", "type": "IBD1", "start_cM": "0", "end_cM": "10", "length": "10"}
    observed = tmp_path / "synthetic_segments.csv"
    with observed.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, row)
        writer.writeheader()
        writer.writerow(row)
    genetic_map = tmp_path / "synthetic.map"
    genetic_map.write_text("1 0 0 0\n1 20 5 15\n")
    records = pair_features([row], 3350.933402, read_sex_map(genetic_map))
    values, features = matrix(records)
    labels = list(IBD_LABEL_MEIOSES)
    model = RandomForestClassifier(n_estimators=2, random_state=1).fit(
        np.repeat(values, len(labels), axis=0), labels,
    )
    release = json.loads((package / RELEASE).read_text())
    release["release_id"] = "synthetic-fixture"
    model_path = tmp_path / "artifacts/task9_releases/synthetic-fixture/model.joblib"
    model_path.parent.mkdir(parents=True)
    joblib.dump({"model": model, "features": features}, model_path)
    release["model"] = {"path": str(model_path.relative_to(tmp_path)), "sha256": digest(model_path)}
    release["features"] = features
    release["sex_map_sha256"] = digest(genetic_map)
    release["environment"]["scikit_learn"] = importlib.metadata.version("scikit-learn")
    write_json(package / RELEASE, release)
    # Construct a separate, self-consistent fixture release for this runtime.
    # No production artifact or evidence is changed; this is not a real release.
    runtime = {
        "python": platform.python_version(), "system": platform.system(), "machine": platform.machine(),
    }
    script = package / SOURCES[0]
    script.write_text(script.read_text().replace(
        'RUNTIME = {"python": "3.11.11", "system": "Linux", "machine": "x86_64"}',
        f"RUNTIME = {runtime!r}",
    ))
    packages = json.loads((package / "package.json").read_text())["runtime"]["packages"]
    (package / REQUIREMENTS).write_text("".join(
        f"{name}=={importlib.metadata.version(name)}\n" for name in packages
    ))
    for name in EVIDENCE[1:]:
        path = package / name
        evidence = json.loads(path.read_text())
        evidence["model_sha256"] = release["model"]["sha256"]
        if name == EVIDENCE[4]:
            for dataset in evidence["datasets"]:
                for case in dataset["cases"].values():
                    case["release_id"] = release["release_id"]
        write_json(path, evidence)
    release["evidence_sha256"] = {name: digest(package / name) for name in EVIDENCE[:4]}
    write_json(package / RELEASE, release)
    fixture_package = tmp_path / "synthetic-package"
    subprocess.run([sys.executable, str(script), "prepare", "--checkout", str(package),
                    "--release", str(package / RELEASE), "--out", str(fixture_package)],
                   capture_output=True, text=True, check=True)
    package = fixture_package
    provenance = tmp_path / "synthetic_provenance.json"
    write_json(provenance, {"entity_id": "syn75202631", "version": 1, "sha256": digest(observed)})
    run_dir = tmp_path / "synthetic-run"
    command = [sys.executable, str(package / "scripts/prepare_task9_submission_package.py"), "run",
               "--package", str(package), "--artifact-root", str(tmp_path), "--sex-map", str(genetic_map),
               "--observed", str(observed), "--input-provenance", str(provenance), "--run-dir", str(run_dir),
               "--expected-package-sha256", digest(package / "package.json")]
    return command, run_dir, observed, provenance, package


def test_synthetic_run_uses_pinned_predictor_and_independent_validator(synthetic_run):
    command, run_dir, _, _, _ = synthetic_run
    result = subprocess.run(command, capture_output=True, text=True, check=True)
    report = json.loads(result.stdout)
    assert report["status"] == "LOCAL_FINAL_FILE_VALIDATED_NOT_SUBMITTED"
    assert report["release_id"] == "synthetic-fixture"
    assert report["validation"]["pairs"] == 1
    assert report["validation"]["input_pair_coverage_verified"] is True
    assert report["submission_count"] == 0
    assert report["output_sha256"] == digest(report["candidate"])
    assert "--release" in json.loads((run_dir / "commands.json").read_text())["predictor"]
    assert "--model" not in json.loads((run_dir / "commands.json").read_text())["predictor"]
    repeated = subprocess.run(command, capture_output=True, text=True)
    assert repeated.returncode != 0
    assert "historical candidates" in repeated.stderr


@pytest.mark.parametrize("change", ["duplicate", "reversed_duplicate", "hbd", "nan", "bad_length", "unknown_chromosome"])
def test_synthetic_invalid_segments_never_produce_a_handoff(synthetic_run, change):
    command, run_dir, observed, provenance, _ = synthetic_run
    with observed.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if change in {"duplicate", "reversed_duplicate"}:
        duplicate = dict(rows[0])
        if change == "reversed_duplicate":
            duplicate["indID1"], duplicate["indID2"] = duplicate["indID2"], duplicate["indID1"]
        rows.append(duplicate)
    else:
        field, value = {"hbd": ("type", "HBD"), "nan": ("end_bp", "nan"),
                        "bad_length": ("length", "11"), "unknown_chromosome": ("chr", "22")}[change]
        rows[0][field] = value
    with observed.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, rows[0])
        writer.writeheader()
        writer.writerows(rows)
    write_json(provenance, {"entity": "syn75202631", "version": 1, "sha256": digest(observed)})
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode != 0
    assert not (run_dir / "handoff.json").exists()
    assert not (run_dir / "prediction/task9_final_candidate.csv").exists()


@pytest.mark.parametrize("mutation", ["model", "map", "roster", "input", "package"])
def test_changes_after_prediction_never_produce_a_handoff(synthetic_run, mutation):
    command, run_dir, _, _, package = synthetic_run
    # A fresh process imports only the fixture package, and intercepts completed
    # subprocess calls to simulate files changing during an otherwise real run.
    harness = '''
import importlib.util, json, pathlib, sys
sys.dont_write_bytecode = True
spec = importlib.util.spec_from_file_location("handoff", sys.argv[1])
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
original = module.subprocess.run
mutation = sys.argv[2]
run_dir = pathlib.Path(sys.argv[3])
package = pathlib.Path(sys.argv[4])
def intercepted(command, **kwargs):
    result = original(command, **kwargs)
    if "validate_submission.py" in command[1]:
        if mutation in {"model", "map"}:
            path = run_dir / "prediction/summary.json"
            summary = module.read_json(path)
            key = next(key for key in summary["source_sha256"]
                       if key.endswith("model.joblib" if mutation == "model" else ".map"))
            summary["source_sha256"][key] = "0" * 64
            module.write_json(path, summary)
        elif mutation == "roster":
            with (run_dir / "expected_pairs.csv").open("a") as handle:
                handle.write("extra_A,extra_B\\n")
        elif mutation == "input":
            with (run_dir / "input.csv").open("a") as handle:
                handle.write("# altered input\\n")
        else:
            with (package / "package.json").open("a") as handle:
                handle.write("\\n")
    return result
module.subprocess.run = intercepted
sys.argv = [sys.argv[1], *sys.argv[5:]]
module.main()
'''
    result = subprocess.run([sys.executable, "-c", harness, command[1], mutation,
                             str(run_dir), str(package), *command[2:]],
                            capture_output=True, text=True)
    assert result.returncode != 0
    assert "validation failed" in result.stderr or "externally pinned checksum" in result.stderr
    assert not (run_dir / "handoff.json").exists()
