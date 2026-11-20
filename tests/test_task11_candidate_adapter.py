"""Boundary checks for the separate truth-free adapter, without detector work."""
import csv
import gzip
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

PATH = Path(__file__).resolve().parents[1] / "scripts/prepare_task11_candidate.py"
SPEC = importlib.util.spec_from_file_location("task11_adapter", PATH)
adapter = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(adapter)


def test_source_rejects_corruption_extra_code_and_escape(tmp_path):
    (tmp_path / "scripts").mkdir()
    path = tmp_path / "scripts/example.py"
    path.write_text("value = 1\n")
    entries = {"scripts/example.py": adapter.digest(path)}
    adapter.verify_source(tmp_path, entries)
    path.write_text("value = 2\n")
    with pytest.raises(ValueError, match="mismatched"):
        adapter.verify_source(tmp_path, entries)
    entries["scripts/example.py"] = adapter.digest(path)
    (tmp_path / "scripts/extra.py").write_text("pass\n")
    with pytest.raises(ValueError, match="extra or missing"):
        adapter.verify_source(tmp_path, entries)
    with pytest.raises(ValueError, match="Unsafe"):
        adapter.verify_source(tmp_path, {"../escape.py": "0" * 64})


def test_source_rejects_symlink(tmp_path):
    (tmp_path / "scripts").mkdir()
    target = tmp_path / "target.py"
    target.write_text("pass\n")
    (tmp_path / "scripts/link.py").symlink_to(target)
    with pytest.raises(ValueError, match="Unsafe"):
        adapter.verify_source(tmp_path, {"scripts/link.py": adapter.digest(target)})


@pytest.mark.parametrize("schema", adapter.SCHEMA_METHODS)
def test_candidate_requires_exact_bound_model_and_source(tmp_path, schema):
    source = tmp_path / "source"
    (source / "scripts").mkdir(parents=True)
    script = source / "scripts/a.py"
    script.write_text("pass\n")
    model = tmp_path / "model.joblib"
    model.write_bytes(b"not loaded by this boundary test")
    protocol = {"schema": schema, "clustering_policy": "average",
                "method_version": adapter.SCHEMA_METHODS[schema], "source_sha256": {"scripts/a.py": adapter.digest(script)},
                "inputs": {}, "gates": {"example": 1}}
    protocol_path = tmp_path / "protocol.json"
    protocol_path.write_text(json.dumps(protocol))
    candidate = {"method_version": protocol["method_version"], "protocol_sha256": adapter.digest(protocol_path),
                 "frozen_before_confirmation": True, "source_sha256": protocol["source_sha256"],
                 "gates": protocol["gates"], "model_path": str(model),
                 "bound_inputs": {str(model): adapter.digest(model)}}
    candidate_path = tmp_path / "candidate.json"
    candidate_path.write_text(json.dumps(candidate))
    arguments = (protocol_path, candidate_path, source, adapter.digest(protocol_path), adapter.digest(candidate_path))
    adapter.verify_candidate(*arguments)
    model.write_bytes(b"corrupt")
    with pytest.raises(ValueError, match="mismatched"):
        adapter.verify_candidate(*arguments)
    model.unlink()
    with pytest.raises(FileNotFoundError):
        adapter.verify_candidate(*arguments)
    candidate["bound_inputs"] = {}
    candidate_path.write_text(json.dumps(candidate))
    with pytest.raises(ValueError, match="not hash-bound"):
        adapter.verify_candidate(protocol_path, candidate_path, source, adapter.digest(protocol_path), adapter.digest(candidate_path))


def write_vcf(path, names=None, body=None):
    names = names or [f"dog,{i}" for i in range(8)]
    body = body or "1\t100\t.\tA\tC\t.\tPASS\t.\tGT:PS\t" + "\t".join(["0|1:7"] * len(names)) + "\n"
    path.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(names) + "\n" + body)
    return body


def test_stage_changes_header_only_and_csv_restores_arbitrary_ids(tmp_path):
    vcf = tmp_path / "input.vcf"
    body = write_vcf(vcf)
    target = tmp_path / "staged.vcf.gz"
    roster = adapter.stage_vcf(vcf, target)
    with gzip.open(target, "rt") as handle:
        assert handle.read().splitlines()[-1] + "\n" == body
    assert roster["internal_samples"] == [f"adapter11_{i}" for i in range(8)]
    rows, families = adapter.restore_predictions(
        [{"indID1": "adapter11_0", "indID2": "adapter11_1", "relation": "siblings"}],
        [["adapter11_0", "adapter11_1"]], roster["internal_to_original"])
    assert families == [["dog,0", "dog,1"]]
    adapter.export_csv(tmp_path / "output.csv", rows)
    with (tmp_path / "output.csv").open() as handle:
        assert list(csv.DictReader(handle)) == rows
    # Compression detection is based on bytes, not a misleading extension.
    again = adapter.stage_vcf(target, tmp_path / "second.gz")
    assert again["variant_records"] == 1
    assert again["unchanged_variant_lines_sha256"] == roster["unchanged_variant_lines_sha256"]


@pytest.mark.parametrize("names", [["same"] * 8, [f"dog {i}" for i in range(8)], [f"dog{i}" for i in range(7)]])
def test_stage_rejects_invalid_roster(tmp_path, names):
    vcf = tmp_path / "input.vcf"
    write_vcf(vcf, names)
    with pytest.raises(ValueError, match="unique whitespace-free"):
        adapter.stage_vcf(vcf, tmp_path / "out.gz")


def test_stage_rejects_changed_column_count(tmp_path):
    vcf = tmp_path / "input.vcf"
    write_vcf(vcf, body="1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|1\n")
    with pytest.raises(ValueError, match="columns"):
        adapter.stage_vcf(vcf, tmp_path / "out.gz")


def test_qualification_cannot_silently_bypass_failed_confirmation(tmp_path):
    args = SimpleNamespace(confirmation_summary=None, confirmation_summary_sha256=None,
                           confirmation_audit=None, confirmation_audit_sha256=None, development_only=False)
    with pytest.raises(ValueError, match="Passing confirmation"):
        adapter.verify_qualification(args, "protocol", "candidate")
    args.development_only = True
    assert adapter.verify_qualification(args, "protocol", "candidate")["development_only"]
    summary = {"protocol_sha256": "protocol", "candidate_sha256": "candidate", "predeclared_criteria_passed": True,
               "checks": dict.fromkeys(adapter.CHECKS, True)}
    audit = {"protocol_sha256": "protocol", "candidate_sha256": "candidate", "all_passed": True,
             "label_and_order_permutations": 160, "unchanged": 160}
    for name, value in (("confirmation_summary", summary), ("confirmation_audit", audit)):
        path = tmp_path / (name + ".json")
        path.write_text(json.dumps(value))
        setattr(args, name, path)
        setattr(args, name + "_sha256", adapter.digest(path))
    args.development_only = False
    assert adapter.verify_qualification(args, "protocol", "candidate")["predeclared_criteria_and_robustness_passed"]
    summary["checks"]["missing_true_pairs"] = False
    args.confirmation_summary.write_text(json.dumps(summary))
    args.confirmation_summary_sha256 = adapter.digest(args.confirmation_summary)
    with pytest.raises(ValueError, match="Passing confirmation"):
        adapter.verify_qualification(args, "protocol", "candidate")


def test_qualification_rejects_changed_hash_and_other_candidate(tmp_path):
    path = tmp_path / "summary.json"
    path.write_text(json.dumps({"protocol_sha256": "other", "candidate_sha256": "candidate"}))
    args = SimpleNamespace(confirmation_summary=path, confirmation_summary_sha256="0" * 64,
                           confirmation_audit=None, confirmation_audit_sha256=None, development_only=True)
    with pytest.raises(ValueError, match="mismatched"):
        adapter.verify_qualification(args, "protocol", "candidate")
    args.confirmation_summary_sha256 = adapter.digest(path)
    with pytest.raises(ValueError, match="identity"):
        adapter.verify_qualification(args, "protocol", "candidate")
