from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import yaml


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inventory", type=Path, default=Path("configs/data_inventory.yaml"))
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    inventory = yaml.safe_load(args.inventory.read_text())
    results = []
    for record in inventory["files"]:
        path = Path(record["path"])
        observed = sha256(path) if path.is_file() else None
        results.append(
            {
                "path": str(path),
                "role": record["role"],
                "exists": path.is_file(),
                "expected_sha256": record["sha256"],
                "observed_sha256": observed,
                "matches": observed == record["sha256"],
            }
        )
    report = {
        "inventory": str(args.inventory),
        "files": len(results),
        "all_match": all(record["matches"] for record in results),
        "results": results,
    }
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    if not report["all_match"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
