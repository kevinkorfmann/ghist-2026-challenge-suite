from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

import yaml


def nested_get(mapping: dict, dotted: str) -> dict:
    value = mapping
    for part in dotted.split("."):
        value = value[part]
    return value


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description="Download registered Synapse inputs")
    parser.add_argument("--challenge", required=True, help="Dotted registry key")
    parser.add_argument("--dataset", choices=("testing", "final"), default="testing")
    parser.add_argument("--registry", default="configs/challenges.yaml")
    parser.add_argument("--version", type=int, help="Pin an entity version; default is the current version")
    parser.add_argument("--data-root", type=Path, default=Path("data/current"))
    args = parser.parse_args()
    token = os.environ.get("SYNAPSE_AUTH_TOKEN")
    if not token:
        raise RuntimeError("SYNAPSE_AUTH_TOKEN is not set")
    import synapseclient
    registry = yaml.safe_load(Path(args.registry).read_text())
    entry = nested_get(registry, args.challenge)
    entity_id = entry[args.dataset]
    syn = synapseclient.Synapse()
    syn.login(authToken=token, silent=True)
    metadata = syn.get(entity_id, version=args.version, downloadFile=False)
    version = int(metadata.versionNumber)
    destination = args.data_root / entry["path"].removeprefix("challenges/") / args.dataset / f"v{version}"
    destination.mkdir(parents=True, exist_ok=True)
    entity = syn.get(entity_id, version=version, downloadLocation=str(destination))
    path = Path(entity.path)
    manifest = {
        "challenge": args.challenge,
        "dataset": args.dataset,
        "entity_id": entity_id,
        "version": version,
        "file_handle_id": str(entity.dataFileHandleId),
        "name": path.name,
        "size_bytes": path.stat().st_size,
        "sha256": sha256(path),
    }
    (destination / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()

