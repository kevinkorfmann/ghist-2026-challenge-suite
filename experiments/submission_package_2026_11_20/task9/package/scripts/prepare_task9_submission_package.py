"""Freeze an offline Task 9 handoff; `run` is a separate, explicit final-data action.

Preparation never opens challenge data or a model. There is no upload operation.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tempfile


SCHEMA = "ghist.task9.submission-package.v1"
SOURCES = (
    "scripts/prepare_task9_submission_package.py",
    "scripts/predict_task9_current.py",
    "scripts/validate_submission.py",
    "src/ghist_relatedness/__init__.py",
    "src/ghist_relatedness/ibd.py",
    "src/ghist_relatedness/labels.py",
    "src/ghist_relatedness/task9_release.py",
    "src/ghist_submission/__init__.py",
    "src/ghist_submission/validate.py",
)
EVIDENCE = (
    "experiments/scientific_2026_11_20/relatedness/task9_current_training.json",
    "experiments/scientific_2026_11_20/relatedness/task9_confirmation_current_interference.json",
    "experiments/scientific_2026_11_20/relatedness/task9_confirmation_current_poisson.json",
    "experiments/readiness_2026_11_20/task9/operational_unversioned.json",
    "experiments/readiness_2026_11_20/task9/operational_versioned.json",
)
REQUIREMENTS = "requirements/task9-release-2026-11-20.txt"
RELEASE = "configs/task9_current_release.json"
RUNTIME = {"python": "3.11.11", "system": "Linux", "machine": "x86_64"}
GATES = {
    "matched_confirmation_at_least_10000_pairs", "both_exact_accuracy_at_least_80pct",
    "both_balanced_accuracy_at_least_75pct", "both_meioses_accuracy_at_least_84pct",
    "every_map_stress_has_probability_support", "exactly_current_training_labels",
    "no_answer_or_confirmation_selection", "distinct_calibration_and_confirmation",
    "same_fixed_model", "production_pair_coverage_and_format",
    "production_matches_fixed_confirmation",
}
LABELS = {
    "avuncular", "double_cousin", "full_1st_cousin", "full_2nd_cousin",
    "full_3rd_cousin", "full_4th_cousin", "full_5th_cousin",
    "grandparent_grandchild", "half_1st_cousin", "siblings",
}


def digest(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text())


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def release_contract(release):
    """Reject historical freezes and default/unversioned model paths."""
    release_id = release.get("release_id", "")
    if (
        not release_id or Path(release_id).name != release_id
        or release_id in {".", ".."}
        or release.get("schema") != "ghist.task9.release.v1"
        or release.get("status") != "VALIDATED_RELEASE_CANDIDATE"
        or set(release.get("accepted_labels", [])) != LABELS
        or len(release.get("accepted_labels", [])) != len(LABELS)
        or release.get("pair_domain") != "unordered_pairs_present_in_segment_csv"
        or set(release.get("gates", {})) != GATES
        or set(release.get("feature_contract_sha256", {})) != {"ibd.py", "labels.py"}
        or not all(value is True for value in release["gates"].values())
        or release.get("model", {}).get("path")
        != f"artifacts/task9_releases/{release_id}/model.joblib"
    ):
        raise ValueError("Require a validated, versioned ten-label release; historical/default paths are refused")


def package_metadata(checkout, release, hashes):
    """Derive metadata from pinned files again during verification."""
    for name, expected in release["feature_contract_sha256"].items():
        if hashes[f"src/ghist_relatedness/{name}"] != expected:
            raise ValueError(f"Release feature source mismatch: {name}")
    training, interference, poisson, operational = (read_json(checkout / name) for name in EVIDENCE[:4])
    # The four release evidence hashes use remote paths; match by content identity.
    if (
        len(release["evidence_sha256"]) != 4
        or set(release["evidence_sha256"].values()) != {hashes[name] for name in EVIDENCE[:4]}
    ):
        raise ValueError("Local evidence differs from the validated release evidence")
    versioned = read_json(checkout / EVIDENCE[4])
    if any(row.get("model_sha256") != release["model"]["sha256"]
           for row in (interference, poisson, operational, versioned)):
        raise ValueError("Confirmation/operational evidence uses a different model")
    if any(row.get("passed") is not True for row in (operational, versioned)):
        raise ValueError("Operational evidence does not pass for this model")
    expected_by_hash = {row["segment_sha256"]: row["runs"]["unchanged"]
                        for row in (interference, poisson)}
    for report in (operational, versioned):
        datasets = report["datasets"]
        if len(datasets) != 2 or {row["source_sha256"] for row in datasets} != set(expected_by_hash):
            raise ValueError("Operational evidence does not cover both confirmations")
        for row in datasets:
            expected = expected_by_hash[row["source_sha256"]]
            if (
                row["expected_pairs"] != expected["pairs"] or row["permutation_invariant"] is not True
                or set(row["cases"]) != {"official_csv", "permuted_reversed_bom_csv"}
                or any(case["validation"]["format_valid"] is not True
                       or case["validation"]["input_pair_coverage_verified"] is not True
                       or case["validation"]["pairs"] != expected["pairs"]
                       or any(not math.isclose(case[metric], expected[metric], rel_tol=0, abs_tol=1e-12)
                              for metric in ("accuracy", "meioses_accuracy"))
                       or (report is versioned and case["release_id"] != release["release_id"])
                       for case in row["cases"].values())
            ):
                raise ValueError("Operational evidence identity/coverage/accuracy mismatch")
    packages = {}
    for line in (checkout / REQUIREMENTS).read_text().splitlines():
        if line.strip() and not line.startswith("#"):
            name, version = line.strip().split("==")
            packages[name] = version
    if packages.get("scikit-learn") != release["environment"]["scikit_learn"]:
        raise ValueError("Release/runtime scikit-learn versions disagree")
    return {
        "schema": SCHEMA,
        "status": "OFFLINE_HANDOFF_FINAL_INPUT_PENDING",
        "release_id": release["release_id"],
        "release_sha256": hashes[RELEASE],
        "files_sha256": hashes,
        "runtime": {**RUNTIME, "packages": packages},
        "final_input": {"entity": "syn75202631", "version": None, "sha256": None, "path": None},
        "final_prediction": {"path": None, "sha256": None, "pairs": None},
        "submission_count": 0,
        "external_artifacts": {"model": release["model"], "sex_map_sha256": release["sex_map_sha256"]},
        "confirmation_pairs": sum(row["runs"]["unchanged"]["pairs"] for row in (interference, poisson)),
        "training_pairs": training["training_pairs"],
        "operational_passed": operational["passed"],
    }


def new_directory(path, description):
    """Check the lexical path before resolve can hide a dangling symlink."""
    path = Path(path).absolute()
    if path.exists() or path.is_symlink():
        raise ValueError(f"{description} already exists; historical candidates must not be reused")
    if any(parent.is_symlink() for parent in path.parents):
        raise ValueError(f"{description} must not have symlink parents")
    return path


def prepare_package(checkout, release_path, out):
    checkout, release_path = map(Path, (checkout, release_path))
    out = new_directory(out, "Package destination")
    release = read_json(release_path)
    release_contract(release)
    sources = {name: checkout / name for name in (*SOURCES, *EVIDENCE, REQUIREMENTS)}
    sources[RELEASE] = release_path
    hashes = {name: digest(path) for name, path in sources.items()}
    package = package_metadata(checkout, release, hashes)
    out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".task9-package-", dir=out.parent) as temporary:
        stage = Path(temporary) / "package"
        stage.mkdir()
        for name, source in sources.items():
            destination = stage / name
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, destination)
            if digest(destination) != hashes[name]:
                raise ValueError(f"Source changed during packaging: {name}")
        write_json(stage / "package.json", package)
        # mkdir is exclusive even if another process creates the destination
        # after the initial check; rename could replace an empty directory.
        out.mkdir()
        for child in stage.iterdir():
            child.rename(out / child.name)
    return package


def verify_package(package_dir, expected_package_sha256=None):
    package_dir = Path(package_dir).resolve()
    if (package_dir / "package.json").is_symlink():
        raise ValueError("Package symlink refused: package.json")
    if expected_package_sha256 is not None and digest(package_dir / "package.json") != expected_package_sha256:
        raise ValueError("Package manifest differs from the externally pinned checksum")
    package = read_json(package_dir / "package.json")
    if package.get("schema") != SCHEMA:
        raise ValueError("Unknown package schema")
    expected_files = set((*SOURCES, *EVIDENCE, REQUIREMENTS, RELEASE))
    if set(package.get("files_sha256", {})) != expected_files:
        raise ValueError("Incomplete package file inventory")
    actual_files = {str(path.relative_to(package_dir)) for path in package_dir.rglob("*")
                    if path.is_file() or path.is_symlink()}
    if actual_files != expected_files | {"package.json"}:
        raise ValueError("Unexpected or missing package files; unpinned bytecode/imports are refused")
    for name, expected in package["files_sha256"].items():
        path = package_dir / name
        if path.is_symlink() or any(parent.is_symlink() for parent in path.parents if parent != package_dir):
            raise ValueError(f"Package symlink refused: {name}")
        if digest(package_dir / name) != expected:
            raise ValueError(f"Package checksum mismatch: {name}")
    release = read_json(package_dir / RELEASE)
    release_contract(release)
    if package != package_metadata(package_dir, release, package["files_sha256"]):
        raise ValueError("Package metadata differs from pinned release/evidence/runtime files")
    return package


def preflight(package_dir, artifact_root, sex_map, expected_package_sha256=None):
    package_dir = Path(package_dir).resolve()
    package = verify_package(package_dir, expected_package_sha256)
    runtime = package["runtime"]
    observed = {"python": platform.python_version(), "system": platform.system(), "machine": platform.machine()}
    if any(observed[name] != runtime[name] for name in observed):
        raise ValueError(f"Use the recorded release runtime: expected {runtime}, observed {observed}")
    for name, version in runtime["packages"].items():
        if importlib.metadata.version(name) != version:
            raise ValueError(f"Release runtime mismatch: {name} must be {version}")
    sys.dont_write_bytecode = True
    sys.path.insert(0, str(package_dir / "src"))
    import joblib
    from ghist_relatedness.task9_release import load_release, validate_payload
    for name, module in tuple(sys.modules.items()):
        if (name == "ghist_relatedness" or name.startswith("ghist_relatedness.")) and not Path(
            module.__file__
        ).resolve().is_relative_to(package_dir / "src"):
            raise ValueError("Preflight must run in a fresh process using the pinned package imports")

    manifest, model_path, _ = load_release(package_dir / RELEASE, Path(artifact_root), Path(sex_map))
    validate_payload(joblib.load(model_path), manifest)
    return {"package_verified": True, "release_id": package["release_id"], "runtime": runtime,
            "model": str(model_path), "model_sha256": manifest["model"]["sha256"],
            "sex_map": str(Path(sex_map).resolve()), "sex_map_sha256": manifest["sex_map_sha256"]}


def run_final(package_dir, artifact_root, sex_map, observed, provenance_path, run_dir, workers=2,
              expected_package_sha256=None):
    """Future action: pin supplied final input, predict once, independently validate.

    No input discovery, download, default candidate reuse, tuning, or submission.
    """
    run_dir = new_directory(run_dir, "Run directory")
    package_dir, artifact_root, sex_map, observed, provenance_path = (
        Path(value).resolve() for value in
        (package_dir, artifact_root, sex_map, observed, provenance_path)
    )
    if workers < 1:
        raise ValueError("workers must be positive")
    provenance = read_json(provenance_path)
    if (
        provenance.get("entity") != "syn75202631"
        or type(provenance.get("version")) is not int or provenance["version"] < 1
        or provenance.get("sha256") != digest(observed)
    ):
        raise ValueError("Final input needs matching syn75202631 entity/version/checksum provenance")
    checked = preflight(package_dir, artifact_root, sex_map, expected_package_sha256)
    initial_package_sha256 = digest(package_dir / "package.json")
    from ghist_relatedness.ibd import read_segment_rows
    from ghist_relatedness.task9_release import input_pairs

    # Copy to a fresh run directory so the predictor and validator use one input.
    run_dir.mkdir(parents=True)
    input_copy = run_dir / "input.csv"
    shutil.copyfile(observed, input_copy)
    if digest(input_copy) != provenance["sha256"]:
        raise ValueError("Input changed during snapshot")
    write_json(run_dir / "input_provenance.json", provenance)
    pairs = input_pairs(read_segment_rows(input_copy))
    with (run_dir / "expected_pairs.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(("indID1", "indID2"))
        writer.writerows(sorted(pairs))
    write_json(run_dir / "preflight.json", checked)
    roster_sha256 = digest(run_dir / "expected_pairs.csv")
    env = dict(os.environ, PYTHONPATH=os.pathsep.join((str(package_dir / "src"), str(package_dir))),
               PYTHONDONTWRITEBYTECODE="1")
    predictor = [sys.executable, str(package_dir / "scripts/predict_task9_current.py"),
                 "--release", str(package_dir / RELEASE), "--artifact-root", str(artifact_root),
                 "--observed", str(input_copy), "--sex-map", str(sex_map), "--dataset", "final",
                 "--input-provenance", str(run_dir / "input_provenance.json"),
                 "--pairs", str(run_dir / "expected_pairs.csv"), "--workers", str(workers),
                 "--out-dir", str(run_dir / "prediction")]
    candidate = run_dir / "prediction/task9_final_candidate.csv"
    validator = [sys.executable, str(package_dir / "scripts/validate_submission.py"), str(candidate),
                 "--challenge", "9", "--dataset", "final", "--pairs", str(run_dir / "expected_pairs.csv"),
                 "--report", str(run_dir / "validation.json")]
    write_json(run_dir / "commands.json", {"predictor": predictor, "validator": validator})
    for name, command in (("predictor", predictor), ("validator", validator)):
        with (run_dir / f"{name}.log").open("w") as log:
            subprocess.run(command, cwd=package_dir, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
    summary = read_json(run_dir / "prediction/summary.json")
    validation = read_json(run_dir / "validation.json")
    if (
        summary["release_id"] != checked["release_id"] or summary["dataset"] != "final"
        or summary["synthetic_input"] is not False or summary["input_provenance"] != provenance
        or summary["source_sha256"][str(input_copy)] != provenance["sha256"]
        or digest(input_copy) != provenance["sha256"]
        or read_json(run_dir / "input_provenance.json") != provenance
        or digest(run_dir / "expected_pairs.csv") != roster_sha256
        or summary["source_sha256"][checked["model"]] != checked["model_sha256"]
        or summary["source_sha256"][str(sex_map)] != checked["sex_map_sha256"]
        or summary["output_sha256"] != digest(candidate)
        or summary["submission_count"] != 0
        or validation["challenge"] != 9 or validation["dataset"] != "final"
        or validation["format_valid"] is not True or validation["input_pair_coverage_verified"] is not True
        or validation["pairs"] != len(pairs)
    ):
        raise ValueError("Final output identity/coverage validation failed")
    verify_package(package_dir, initial_package_sha256)
    report = {"status": "LOCAL_FINAL_FILE_VALIDATED_NOT_SUBMITTED", "submission_count": 0,
              "release_id": checked["release_id"], "input_provenance": provenance,
              "package_sha256": digest(package_dir / "package.json"), "validation": validation,
              "candidate": str(candidate), "output_sha256": digest(candidate),
              "run_files_sha256": {str(path.relative_to(run_dir)): digest(path)
                                   for path in sorted(run_dir.rglob("*")) if path.is_file()}}
    write_json(run_dir / "handoff.json", report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    actions = parser.add_subparsers(dest="action", required=True)
    prepare = actions.add_parser("prepare", help="Offline: snapshot source/evidence only")
    prepare.add_argument("--checkout", type=Path, required=True)
    prepare.add_argument("--release", type=Path, required=True)
    prepare.add_argument("--out", type=Path, required=True)
    verify = actions.add_parser("verify", help="Offline: verify snapshot checksums only")
    verify.add_argument("--package", type=Path, required=True)
    verify.add_argument("--expected-package-sha256")
    for action in ("preflight", "run"):
        command = actions.add_parser(action, help="Verify model/runtime" if action == "preflight" else "Run supplied final input locally; never submit")
        command.add_argument("--package", type=Path, required=True)
        command.add_argument("--artifact-root", type=Path, required=True)
        command.add_argument("--sex-map", type=Path, required=True)
        command.add_argument("--expected-package-sha256", required=True,
                             help="Manifest SHA-256 from the trusted handoff, not the transferred manifest")
        if action == "run":
            command.add_argument("--observed", type=Path, required=True)
            command.add_argument("--input-provenance", type=Path, required=True)
            command.add_argument("--run-dir", type=Path, required=True)
            command.add_argument("--workers", type=int, default=2)
    args = parser.parse_args()
    if args.action == "prepare":
        result = prepare_package(args.checkout, args.release, args.out)
    elif args.action == "verify":
        result = verify_package(args.package, args.expected_package_sha256)
    elif args.action == "preflight":
        result = preflight(args.package, args.artifact_root, args.sex_map, args.expected_package_sha256)
    else:
        result = run_final(args.package, args.artifact_root, args.sex_map, args.observed,
                           args.input_provenance, args.run_dir, args.workers, args.expected_package_sha256)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
