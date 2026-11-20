import hashlib
import json

import pytest
import yaml

from scripts import verify_data_inventory


def make_inventory(tmp_path, records):
    inventory = tmp_path / "inventory.yaml"
    inventory.write_text(yaml.safe_dump({"files": records}))
    return inventory


def file_record(tmp_path, name="input", content=b"fixture"):
    path = tmp_path / name
    path.write_bytes(content)
    return {"path": str(path), "role": "testing", "sha256": hashlib.sha256(content).hexdigest()}


def test_valid_inventory_checks_bytes_and_accepts_uppercase_digest(tmp_path):
    record = file_record(tmp_path)
    record["sha256"] = record["sha256"].upper()
    report = verify_data_inventory.verify_inventory(make_inventory(tmp_path, [record]))
    assert report["all_match"]
    assert report["results"][0]["observed_sha256"] == record["sha256"].lower()


@pytest.mark.parametrize("expected", [None, "", "not-a-hash", "0" * 63, 123])
def test_missing_file_and_invalid_hash_never_match(tmp_path, expected):
    inventory = make_inventory(tmp_path, [{"path": str(tmp_path / "missing"), "sha256": expected}])
    report = verify_data_inventory.verify_inventory(inventory)
    assert not report["all_match"]
    assert not report["results"][0]["matches"]
    assert "hex digest" in report["results"][0]["error"]


@pytest.mark.parametrize("contents", ["files: []\n", "{}\n", "null\n", "files: nope\n", "[invalid yaml"])
def test_empty_or_malformed_inventory_fails_closed(tmp_path, contents):
    inventory = tmp_path / "inventory.yaml"
    inventory.write_text(contents)
    report = verify_data_inventory.verify_inventory(inventory)
    assert not report["all_match"]
    assert report["errors"]


def test_read_failure_and_missing_files_do_not_hide_later_results(tmp_path, monkeypatch):
    blocked = file_record(tmp_path, "blocked")
    valid = file_record(tmp_path, "valid")
    missing = {**valid, "path": str(tmp_path / "missing")}
    directory = {**valid, "path": str(tmp_path)}
    mismatch = {**valid, "sha256": "0" * 64}
    original_sha256 = verify_data_inventory.sha256

    def read(path):
        if path.name == "blocked":
            raise PermissionError("unreadable fixture")
        return original_sha256(path)

    monkeypatch.setattr(verify_data_inventory, "sha256", read)
    report = verify_data_inventory.verify_inventory(
        make_inventory(tmp_path, [blocked, missing, directory, mismatch, "bad entry", {}, valid])
    )
    assert not report["all_match"]
    assert report["files"] == 7
    assert [row["matches"] for row in report["results"]] == [False] * 6 + [True]
    assert "PermissionError" in report["results"][0]["error"]
    assert "missing" in report["results"][1]["error"]
    assert "SHA-256 mismatch" == report["results"][3]["error"]


def test_cli_writes_complete_failure_report_and_exits_nonzero(tmp_path, monkeypatch, capsys):
    record = file_record(tmp_path)
    inventory = make_inventory(tmp_path, [record, {**record, "path": str(tmp_path / "missing")}])
    output = tmp_path / "reports" / "audit.json"
    monkeypatch.setattr("sys.argv", ["verify_data_inventory", "--inventory", str(inventory), "--report", str(output)])
    with pytest.raises(SystemExit) as raised:
        verify_data_inventory.main()
    assert raised.value.code == 1
    report = json.loads(output.read_text())
    assert report == json.loads(capsys.readouterr().out)
    assert report["files"] == 2
    assert report["results"][0]["matches"]
    assert not report["results"][1]["matches"]


def test_malformed_yaml_fields_and_nul_path_remain_reportable(tmp_path, monkeypatch, capsys):
    inventory = tmp_path / "inventory.yaml"
    inventory.write_text('files:\n  - {path: "bad\\0path", sha256: null}\n  - {path: [bad], sha256: 2026-11-20}\n')
    output = tmp_path / "report.json"
    monkeypatch.setattr("sys.argv", ["verify_data_inventory", "--inventory", str(inventory), "--report", str(output)])
    with pytest.raises(SystemExit) as raised:
        verify_data_inventory.main()
    assert raised.value.code == 1
    report = json.loads(output.read_text())
    assert report == json.loads(capsys.readouterr().out)
    assert not report["all_match"]
    assert all(row["error"] for row in report["results"])


@pytest.mark.parametrize("target", ["inventory", "input", "input_symlink"])
def test_report_cannot_overwrite_inventory_or_input(tmp_path, target):
    record = file_record(tmp_path)
    inventory = make_inventory(tmp_path, [record])
    report = verify_data_inventory.verify_inventory(inventory)
    destination = inventory if target == "inventory" else tmp_path / "input"
    if target == "input_symlink":
        destination = tmp_path / "link"
        destination.symlink_to(tmp_path / "input")
    original_inventory = inventory.read_bytes()
    with pytest.raises(ValueError, match="must not overwrite"):
        verify_data_inventory.write_report(report, destination)
    assert inventory.read_bytes() == original_inventory
    assert (tmp_path / "input").read_bytes() == b"fixture"


def test_failed_report_replacement_preserves_previous_report(tmp_path, monkeypatch):
    report = verify_data_inventory.verify_inventory(make_inventory(tmp_path, [file_record(tmp_path)]))
    output = tmp_path / "report.json"
    output.write_bytes(b"previous report")
    original_files = set(tmp_path.iterdir())

    def fail_replace(*args):
        raise OSError("replacement failed")

    monkeypatch.setattr(verify_data_inventory.os, "replace", fail_replace)
    with pytest.raises(OSError, match="replacement failed"):
        verify_data_inventory.write_report(report, output)
    assert output.read_bytes() == b"previous report"
    assert set(tmp_path.iterdir()) == original_files
