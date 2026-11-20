from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import tempfile

import yaml


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_record(record) -> dict:
    result = {
        "path": None,
        "role": None,
        "exists": None,
        "expected_sha256": None,
        "observed_sha256": None,
        "matches": False,
    }
    if not isinstance(record, dict):
        result["error"] = "Inventory entry must be a mapping"
        return result
    # YAML scalars such as dates, and malformed nested fields, must still be reportable as JSON.
    for field, key in [("path", "path"), ("role", "role"), ("expected_sha256", "sha256")]:
        value = record.get(key)
        result[field] = value if value is None or isinstance(value, str) else repr(value)
    path_text = record.get("path")
    if not isinstance(path_text, str) or not path_text or "\0" in path_text:
        result["error"] = "Inventory entry must have a nonempty path without NUL characters"
        return result
    path = Path(path_text)
    expected = record.get("sha256")
    if not isinstance(expected, str) or not re.fullmatch(r"[0-9a-fA-F]{64}", expected):
        result["error"] = "Expected SHA-256 must be a 64-character hex digest"
        return result
    try:
        result["exists"] = path.is_file()
        if not result["exists"]:
            result["error"] = "File is missing or is not a regular file"
        else:
            result["observed_sha256"] = sha256(path)
            result["matches"] = result["observed_sha256"] == expected.lower()
            if not result["matches"]:
                result["error"] = "SHA-256 mismatch"
    except (OSError, ValueError) as error:
        result["error"] = f"{type(error).__name__}: {error}"
    return result


def verify_inventory(inventory_path: Path) -> dict:
    errors = []
    records = []
    try:
        inventory = yaml.safe_load(inventory_path.read_text())
        if not isinstance(inventory, dict) or not isinstance(inventory.get("files"), list):
            errors.append("Inventory must contain a files list")
        else:
            records = inventory["files"]
            if not records:
                errors.append("Inventory files list is empty; no files were verified")
    except (OSError, ValueError, yaml.YAMLError) as error:
        errors.append(f"{type(error).__name__}: {error}")
    results = [verify_record(record) for record in records]
    return {
        "inventory": str(inventory_path),
        "files": len(results),
        "all_match": bool(results) and not errors and all(record["matches"] for record in results),
        "errors": errors,
        "results": results,
    }


def write_report(report: dict, destination: Path) -> None:
    protected = {Path(report["inventory"]).resolve()}
    for row in report["results"]:
        if isinstance(row["path"], str):
            try:
                protected.add(Path(row["path"]).resolve())
            except (OSError, RuntimeError, ValueError):
                # Invalid input paths have already failed verification; report those errors.
                pass
    if destination.resolve() in protected:
        raise ValueError("Report path must not overwrite the inventory or an inventoried input")
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", dir=destination.parent, delete=False) as stream:
            temporary = Path(stream.name)
            stream.write(json.dumps(report, indent=2) + "\n")
        os.replace(temporary, destination)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inventory", type=Path, default=Path("configs/data_inventory.yaml"))
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    report = verify_inventory(args.inventory)
    if args.report:
        write_report(report, args.report)
    print(json.dumps(report, indent=2))
    if not report["all_match"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
