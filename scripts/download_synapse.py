from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

import synapseclient
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
    args = parser.parse_args()
    token = os.environ.get("SYNAPSE_AUTH_TOKEN")
    if not token:
        raise RuntimeError("SYNAPSE_AUTH_TOKEN is not set")
    registry = yaml.safe_load(Path(args.registry).read_text())
    entry = nested_get(registry, args.challenge)
    entity_id = entry[args.dataset]
    destination = Path("data/raw") / entry["path"].removeprefix("challenges/") / args.dataset
    destination.mkdir(parents=True, exist_ok=True)
    syn = synapseclient.Synapse()
    syn.login(authToken=token, silent=True)
    entity = syn.get(entity_id, downloadLocation=str(destination))
    path = Path(entity.path)
    manifest = {
        "challenge": args.challenge,
        "dataset": args.dataset,
        "entity_id": entity_id,
        "name": path.name,
        "size_bytes": path.stat().st_size,
        "sha256": sha256(path),
    }
    (destination / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()

