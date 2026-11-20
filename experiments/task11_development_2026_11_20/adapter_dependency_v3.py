"""Truth-free Task 11 CSV export using an immutable, hash-bound candidate.

This adapter does not select a model or change the frozen scientific source.
Without passing confirmation and robustness evidence, --development-only is
required. It never uploads anything or authorizes final challenge access.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import importlib
import json
import os
import subprocess
import sys
from pathlib import Path

VERSION = "task11-truth-free-adapter-v3"
SCHEMA_METHODS = {
    "task11-prospective-development-v2": "task11-detected-ibd-pair-first-development-v2",
    "task11-prospective-development-v3": "task11-detected-ibd-pair-first-development-v3",
    "task11-prospective-development-v4": "task11-detected-ibd-pair-first-development-v4",
}
CHECKS = {"accuracy", "balanced_accuracy", "each_label_recall", "bootstrap_lower",
          "family_recovery", "missing_samples", "missing_true_pairs", "extra_pairs", "calibration"}


def digest(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            value.update(block)
    return value.hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def require_hash(path, expected):
    if not isinstance(expected, str) or len(expected) != 64 or digest(path) != expected:
        raise ValueError(f"Missing or mismatched bound file: {path}")


def verify_source(source, entries):
    """Check before importing any candidate code; disallow path escapes/symlinks."""
    source = Path(source).resolve(strict=True)
    if not isinstance(entries, dict) or not entries:
        raise ValueError("Empty or invalid source binding")
    for name, expected in entries.items():
        relative = Path(name)
        path = source / relative
        if (relative.is_absolute() or ".." in relative.parts or
                not path.resolve().is_relative_to(source) or path.is_symlink()):
            raise ValueError(f"Unsafe source path: {name}")
        require_hash(path, expected)
    found = {str(p.relative_to(source)) for p in (source / "scripts").glob("*.py")}
    found.update(str(p.relative_to(source)) for p in (source / "src").rglob("*.py"))
    if found != set(entries):
        raise ValueError("Snapshot has extra or missing Python source files")


def verify_candidate(protocol_path, candidate_path, source, protocol_hash, candidate_hash):
    require_hash(protocol_path, protocol_hash)
    require_hash(candidate_path, candidate_hash)
    protocol, candidate = load(protocol_path), load(candidate_path)
    if (protocol.get("schema") not in SCHEMA_METHODS or
            protocol.get("method_version") != SCHEMA_METHODS[protocol["schema"]] or
            protocol.get("clustering_policy") != "average"):
        raise ValueError("Unsupported protocol or clustering method")
    if (candidate.get("method_version") != protocol.get("method_version") or
            candidate.get("protocol_sha256") != protocol_hash or
            candidate.get("frozen_before_confirmation") is not True or
            candidate.get("source_sha256") != protocol.get("source_sha256") or
            candidate.get("gates") != protocol.get("gates")):
        raise ValueError("Candidate/protocol identity mismatch")
    verify_source(source, protocol["source_sha256"])
    model = str(Path(candidate["model_path"]).resolve())
    if model not in candidate["bound_inputs"]:
        raise ValueError("Model is not hash-bound")
    for path, expected in {**protocol["inputs"], **candidate["bound_inputs"]}.items():
        require_hash(path, expected)
    return protocol, candidate


def verify_qualification(args, protocol_hash, candidate_hash):
    """Qualification evidence is bound independently of this inference input."""
    evidence = {}
    for name in ("confirmation_summary", "confirmation_audit"):
        path, expected = getattr(args, name), getattr(args, name + "_sha256")
        if path is not None:
            require_hash(path, expected)
            value = load(path)
            if value.get("protocol_sha256") != protocol_hash or value.get("candidate_sha256") != candidate_hash:
                raise ValueError(f"Qualification identity mismatch: {name}")
            evidence[name] = value
    summary, audit = evidence.get("confirmation_summary", {}), evidence.get("confirmation_audit", {})
    checks = summary.get("checks", {})
    passed = (summary.get("predeclared_criteria_passed") is True and
              set(checks) == CHECKS and all(value is True for value in checks.values()) and
              audit.get("all_passed") is True and
              audit.get("label_and_order_permutations") == 160 and audit.get("unchanged") == 160)
    if not passed and not args.development_only:
        raise ValueError("Passing confirmation and robustness evidence required; use --development-only for synthetic development")
    return {"predeclared_criteria_and_robustness_passed": passed,
            "development_only": args.development_only or not passed,
            "scope": "Prespecified synthetic qualification; not proof of biological transfer or permission to access final challenge data"}


def stage_vcf(source, target):
    """Rename the roster only, preserving all variant/FORMAT/phase data verbatim."""
    source, target = Path(source), Path(target)
    with source.open("rb") as handle:
        compressed = handle.read(2) == b"\x1f\x8b"
    opener = gzip.open if compressed else open
    samples, data_hash, records = None, hashlib.sha256(), 0
    with opener(source, "rt", encoding="utf-8", newline="") as reader, gzip.open(target, "wt", encoding="utf-8", newline="") as writer:
        for line in reader:
            if line.startswith("#CHROM\t"):
                header = line.rstrip("\r\n").split("\t")
                if samples is not None or header[:9] != ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]:
                    raise ValueError("Invalid or duplicate VCF roster")
                samples = header[9:]
                if len(samples) < 8 or len(set(samples)) != len(samples) or any(not name or any(c.isspace() for c in name) for name in samples):
                    raise ValueError("Expected at least eight unique whitespace-free VCF sample IDs")
                internal = [f"adapter11_{i}" for i in range(len(samples))]
                writer.write("\t".join(header[:9] + internal) + "\n")
            elif line.startswith("##") and samples is None:
                writer.write(line)
            else:
                if samples is None or line.startswith("#"):
                    raise ValueError("Missing roster or misplaced VCF header")
                if len(line.rstrip("\r\n").split("\t")) != 9 + len(samples):
                    raise ValueError("VCF genotype columns differ from the roster")
                data_hash.update(line.encode("utf-8"))
                records += 1
                writer.write(line)
    if samples is None or not records:
        raise ValueError("VCF has no roster or variant records")
    return {"samples": samples, "internal_samples": internal,
            "internal_to_original": dict(zip(internal, samples, strict=True)),
            "variant_records": records, "unchanged_variant_lines_sha256": data_hash.hexdigest(),
            "scope": "Only sample-header aliases change; strict frozen preparation then checks phase and drops incomplete sites"}


def restore_predictions(rows, families, mapping):
    restored = []
    for row in rows:
        if row["indID1"] not in mapping or row["indID2"] not in mapping:
            raise ValueError("Prediction contains an unknown internal sample")
        restored.append({**row, "indID1": mapping[row["indID1"]], "indID2": mapping[row["indID2"]]})
    return restored, [[mapping[s] for s in family] for family in families]


def export_csv(path, rows):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["indID1", "indID2", "relation"])
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row[name] for name in writer.fieldnames})


def run(args):
    protocol, candidate = verify_candidate(args.protocol, args.candidate, args.source,
                                          args.protocol_sha256, args.candidate_sha256)
    qualification = verify_qualification(args, args.protocol_sha256, args.candidate_sha256)
    if not os.environ.get("SLURM_JOB_ID"):
        raise ValueError("Heavy detector inference requires a Slurm allocation")
    for path in (args.map, args.hapibd):
        if str(path.resolve()) not in protocol["inputs"]:
            raise ValueError(f"Unregistered map or detector: {path}")
    source = args.source.resolve()
    # CLI enters Python isolated mode first. Candidate subprocesses inherit only
    # this source tree via PYTHONPATH and run with cwd equal to the snapshot.
    if any(name == "scripts" or name.startswith(("scripts.", "ghist_")) for name in sys.modules):
        raise ValueError("Candidate modules were imported before source verification")
    sys.dont_write_bytecode = True
    paths = [str(source), str(source / "src"), str(source / "scripts")]
    sys.path[:0] = paths
    os.environ["PYTHONPATH"] = os.pathsep.join(paths)
    frozen = importlib.import_module("scripts.develop_task11_pedigree")
    validator = importlib.import_module("ghist_submission.validate")
    frozen.validate_protocol(protocol)
    frozen.require_runtime(protocol)
    input_hash = digest(args.vcf)
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    audit = {"adapter_version": VERSION, "adapter_sha256": digest(__file__),
             "protocol_sha256": args.protocol_sha256, "candidate_sha256": args.candidate_sha256,
             "input_vcf": str(args.vcf.resolve()), "input_vcf_sha256": input_hash,
             "qualification": qualification, "source": str(source), "truth_used": False,
             "job_id": os.environ["SLURM_JOB_ID"], "submission_count": 0,
             "submission_ready": False, "completed": False}
    try:
        # The historical internal filename is fixed by the immutable helper;
        # it is just a staged input alias, with no assertion of synthetic origin.
        roster = stage_vcf(args.vcf, out / "synthetic.vcf.gz")
        dump(out / "roster.json", roster)
        inference = frozen.detect_cluster_extract(out, args.map.resolve(), args.hapibd.resolve(),
                                                  roster["internal_samples"], "average")
        features = frozen.np.load(out / "features.npz")["features"]
        model = frozen.joblib.load(candidate["model_path"])
        if set(model.classes_) != set(protocol["labels"]) or model.n_features_in_ != features.shape[1]:
            raise ValueError("Model labels or feature dimension differs from candidate")
        rows = frozen.predict(load(out / "feature_pairs.json"), features, model,
                              candidate["temperature"], candidate["sibling_pooling"])
        rows, families = restore_predictions(rows, inference["clusters"]["families"], roster["internal_to_original"])
        expected = validator.within_family_pairs(families, roster["samples"])
        export_csv(out / "task11_candidate.csv", rows)
        validation = validator.validate_submission(out / "task11_candidate.csv", 11, args.dataset, expected)
        dump(out / "predictions.json", rows)
        dump(out / "families.json", {"families": families, "supplied_family_count": 8})
        dump(out / "validation.json", validation)
        require_hash(args.vcf, input_hash)
        verify_candidate(args.protocol, args.candidate, source, args.protocol_sha256, args.candidate_sha256)
        verify_qualification(args, args.protocol_sha256, args.candidate_sha256)
        frozen.require_runtime(protocol)
        audit.update(completed=True, samples=len(roster["samples"]), pairs=len(rows),
                     validation=validation, imported_sources=frozen.imported_source_tree(),
                     model_sha256=digest(candidate["model_path"]), stages=inference["stages"],
                     confirmation_evidence_sha256={name: digest(getattr(args, name)) for name in
                         ("confirmation_summary", "confirmation_audit") if getattr(args, name)})
    except Exception as error:
        audit["failure"] = {"type": type(error).__name__, "message": str(error)}
        raise
    finally:
        audit["output_sha256"] = {str(path.relative_to(out)): digest(path) for path in sorted(out.rglob("*")) if path.is_file() and path.name != "adapter_audit.json"}
        dump(out / "adapter_audit.json", audit)
    return audit


def main():
    if not sys.flags.isolated:
        return subprocess.call([sys.executable, "-I", str(Path(__file__).resolve()), *sys.argv[1:]])
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("vcf", "source", "protocol", "candidate", "map", "hapibd", "out"):
        parser.add_argument("--" + name, type=Path, required=True)
    for name in ("protocol", "candidate"):
        parser.add_argument("--" + name + "-sha256", required=True)
    for name in ("confirmation-summary", "confirmation-audit"):
        parser.add_argument("--" + name, type=Path)
        parser.add_argument("--" + name + "-sha256")
    parser.add_argument("--development-only", action="store_true")
    parser.add_argument("--dataset", choices=("testing", "final"), default="testing",
                        help="Validator format context only; does not authorize access to any data")
    args = parser.parse_args()
    if args.development_only and args.dataset != "testing":
        parser.error("Development-only exports require --dataset testing")
    audit = run(args)
    print(json.dumps({"completed": audit["completed"], "pairs": audit["pairs"], "qualification": audit["qualification"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
