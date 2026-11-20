from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import tempfile

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


def sha256_argument(value: str) -> str:
    if not re.fullmatch(r"[0-9a-fA-F]{64}", value):
        raise argparse.ArgumentTypeError("expected a 64-character SHA-256 hex digest")
    return value.lower()


def check_identity(entity, entity_id: str, version: int, file_handle_id: str) -> None:
    if (
        str(getattr(entity, "id", "")) != entity_id
        or getattr(entity, "versionNumber", None) != version
        or str(getattr(entity, "dataFileHandleId", "")) != file_handle_id
    ):
        raise RuntimeError("Downloaded entity identity differs from the resolved version")


def download_version(syn, entity_id: str, entry: dict, args) -> dict:
    relative_path = Path(entry["path"])
    if (
        relative_path.is_absolute()
        or ".." in relative_path.parts
        or len(relative_path.parts) < 2
        or relative_path.parts[0] != "challenges"
    ):
        raise ValueError("Registry path must be a relative path beneath challenges/")
    metadata = syn.get(entity_id, version=args.version, downloadFile=False)
    version = getattr(metadata, "versionNumber", None)
    file_handle_id = str(getattr(metadata, "dataFileHandleId", ""))
    if type(version) is not int or version < 1 or not file_handle_id.isdigit():
        raise RuntimeError("Resolved entity is missing a valid version or file handle")
    check_identity(metadata, entity_id, version, file_handle_id)
    if args.version is not None and args.version != version:
        raise RuntimeError("Resolved entity does not match the requested version")

    destination = args.data_root / Path(*relative_path.parts[1:]) / args.dataset / f"v{version}"
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(f"Refusing to replace existing download: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=f".v{version}-", dir=destination.parent) as temporary:
        staging = Path(temporary).resolve()
        entity = syn.get(entity_id, version=version, downloadLocation=str(staging))
        check_identity(entity, entity_id, version, file_handle_id)
        returned_path = getattr(entity, "path", None)
        if not returned_path:
            raise RuntimeError("Synapse did not return a downloaded file")
        path = Path(returned_path)
        if path.is_symlink() or path.resolve().parent != staging or not path.is_file():
            raise RuntimeError("Synapse did not return a regular file in the staging directory")
        if path.name == "manifest.json":
            raise RuntimeError("Downloaded filename conflicts with manifest.json")

        size = path.stat().st_size
        expected_size = getattr(metadata, "fileSize", None)
        if expected_size is not None and size != int(expected_size):
            raise RuntimeError("Downloaded size differs from Synapse metadata")
        expected_md5 = getattr(metadata, "md5", None)
        if expected_md5:
            with path.open("rb") as stream:
                observed_md5 = hashlib.file_digest(stream, "md5").hexdigest()
            if observed_md5 != str(expected_md5).lower():
                raise RuntimeError("Downloaded MD5 differs from Synapse metadata")
        digest = sha256(path)
        if args.expected_sha256 is not None and digest != args.expected_sha256:
            raise RuntimeError("Downloaded SHA-256 differs from --expected-sha256")
        manifest = {
            "challenge": args.challenge,
            "dataset": args.dataset,
            "entity_id": entity_id,
            "version": version,
            "file_handle_id": file_handle_id,
            "name": path.name,
            "size_bytes": size,
            "sha256": digest,
        }
        if expected_md5:
            manifest["synapse_md5"] = str(expected_md5).lower()
        if args.expected_sha256 is not None:
            manifest["expected_sha256"] = args.expected_sha256
        manifest_path = staging / "manifest.json"
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")

        # Exclusive creation also protects a competing downloader's publication.
        # Publish the completion manifest last; incomplete directories are never reused.
        destination.mkdir(exist_ok=False)
        try:
            path.rename(destination / path.name)
            manifest_path.rename(destination / "manifest.json")
        except BaseException:
            shutil.rmtree(destination)
            raise
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description="Download registered Synapse inputs")
    parser.add_argument("--challenge", required=True, help="Dotted registry key")
    parser.add_argument("--dataset", choices=("testing", "final"), default="testing")
    parser.add_argument("--registry", default="configs/challenges.yaml")
    parser.add_argument("--version", type=int, help="Pin an entity version; default is the current version")
    parser.add_argument("--data-root", type=Path, default=Path("data/current"))
    parser.add_argument("--expected-sha256", type=sha256_argument, help="Verify a separately trusted digest")
    args = parser.parse_args()
    if args.version is not None and args.version < 1:
        parser.error("--version must be positive")
    token = os.environ.get("SYNAPSE_AUTH_TOKEN")
    if not token:
        raise RuntimeError("SYNAPSE_AUTH_TOKEN is not set")
    import synapseclient
    registry = yaml.safe_load(Path(args.registry).read_text())
    entry = nested_get(registry, args.challenge)
    entity_id = entry[args.dataset]
    syn = synapseclient.Synapse()
    syn.login(authToken=token, silent=True)
    manifest = download_version(syn, entity_id, entry, args)
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()

