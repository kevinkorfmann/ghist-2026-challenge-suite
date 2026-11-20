"""Byte-level bindings for independent relatedness validation records."""

from __future__ import annotations

import hashlib
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_hashes(expected: dict[str, str | None], *, context: str) -> None:
    """Reject missing or changed files; never replace an expected digest."""
    for filename, digest in expected.items():
        if (
            not isinstance(digest, str)
            or len(digest) != 64
            or any(char not in "0123456789abcdef" for char in digest)
        ):
            raise ValueError(f"{context}: missing valid recorded SHA-256 for {filename}")
        path = Path(filename)
        if not path.is_file() or sha256(path) != digest:
            raise ValueError(f"{context}: SHA-256 mismatch for {filename}")


def bind_outputs(directory: Path, result: dict, *, protocol_name: str) -> dict:
    """Bind a new result to the protocol and every retained generated file.

    This attests bytes present at finalization; it cannot retrospectively prove
    the identity of files from historical runs that did not record digests.
    """
    return {
        **result,
        "provenance_schema": "relatedness-pipeline-v2",
        "protocol_sha256": sha256(directory / protocol_name),
        "output_hashes": {
            str(path.relative_to(directory)): sha256(path)
            for path in sorted(directory.rglob("*"))
            if path.is_file() and path != directory / "result.json"
        },
    }
