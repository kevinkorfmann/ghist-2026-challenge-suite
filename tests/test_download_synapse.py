import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from scripts import download_synapse


def test_downloader_pins_version_and_keeps_archived_inputs_separate(
    tmp_path, monkeypatch
):
    calls = []

    class Synapse:
        def login(self, **kwargs):
            pass

        def get(self, entity_id, **kwargs):
            calls.append(kwargs)
            if kwargs.get("downloadFile") is False:
                return SimpleNamespace(id=entity_id, versionNumber=2, dataFileHandleId="123")
            destination = tmp_path / kwargs["downloadLocation"] / "input.vcf.gz"
            destination.write_bytes(b"fixture")
            return SimpleNamespace(
                id=entity_id, versionNumber=2, path=str(destination), dataFileHandleId="123"
            )

    registry = tmp_path / "registry.yaml"
    registry.write_text(
        "demography:\n  bottleneck:\n    path: challenges/demography/bottleneck\n    testing: syn123\n"
    )
    monkeypatch.setitem(
        __import__("sys").modules, "synapseclient", SimpleNamespace(Synapse=Synapse)
    )
    monkeypatch.setenv("SYNAPSE_AUTH_TOKEN", "fixture-not-a-real-token")
    monkeypatch.setattr(
        "sys.argv",
        [
            "download_synapse",
            "--challenge",
            "demography.bottleneck",
            "--registry",
            str(registry),
            "--data-root",
            str(tmp_path / "current"),
        ],
    )
    download_synapse.main()
    assert calls[1]["version"] == 2
    manifest = json.loads(
        (
            tmp_path / "current/demography/bottleneck/testing/v2/manifest.json"
        ).read_text()
    )
    assert manifest["version"] == 2
    assert manifest["entity_id"] == "syn123"
    assert manifest["file_handle_id"] == "123"


@pytest.fixture
def download_setup(tmp_path, monkeypatch):
    state = SimpleNamespace(
        calls=[],
        metadata=SimpleNamespace(
            id="syn123", versionNumber=2, dataFileHandleId="123",
            fileSize=7, md5=hashlib.md5(b"fixture").hexdigest(),
        ),
        entity_updates={},
        name="input.vcf.gz",
    )

    def download(staging):
        path = staging / state.name
        path.write_bytes(b"fixture")
        return path

    state.download = download

    class Synapse:
        def login(self, **kwargs):
            pass

        def get(self, entity_id, **kwargs):
            state.calls.append(kwargs)
            if kwargs.get("downloadFile") is False:
                return state.metadata
            path = state.download(Path(kwargs["downloadLocation"]))
            fields = dict(id=entity_id, versionNumber=2, dataFileHandleId="123", path=path)
            fields.update(state.entity_updates)
            return SimpleNamespace(**fields)

    registry = tmp_path / "registry.yaml"
    registry.write_text(
        "demography:\n  bottleneck:\n    path: challenges/demography/bottleneck\n    testing: syn123\n"
    )
    monkeypatch.setitem(__import__("sys").modules, "synapseclient", SimpleNamespace(Synapse=Synapse))
    monkeypatch.setenv("SYNAPSE_AUTH_TOKEN", "fixture-not-a-real-token")
    state.argv = [
        "download_synapse", "--challenge", "demography.bottleneck",
        "--registry", str(registry), "--data-root", str(tmp_path / "current"),
    ]
    monkeypatch.setattr("sys.argv", state.argv)
    state.destination = tmp_path / "current/demography/bottleneck/testing/v2"
    state.registry = registry
    return state


def test_download_verifies_and_records_trusted_hash(download_setup):
    state = download_setup
    digest = hashlib.sha256(b"fixture").hexdigest()
    state.argv.extend(["--expected-sha256", digest.upper()])
    download_synapse.main()
    manifest = json.loads((state.destination / "manifest.json").read_text())
    assert manifest["sha256"] == manifest["expected_sha256"] == digest
    assert manifest["synapse_md5"] == state.metadata.md5
    assert (state.destination / manifest["name"]).read_bytes() == b"fixture"
    assert list(state.destination.parent.iterdir()) == [state.destination]


def test_existing_download_and_manifest_are_preserved(download_setup):
    state = download_setup
    state.destination.mkdir(parents=True)
    (state.destination / "input.vcf.gz").write_bytes(b"archived bytes")
    (state.destination / "manifest.json").write_bytes(b"archived manifest")
    with pytest.raises(FileExistsError, match="Refusing to replace"):
        download_synapse.main()
    assert len(state.calls) == 1
    assert (state.destination / "input.vcf.gz").read_bytes() == b"archived bytes"
    assert (state.destination / "manifest.json").read_bytes() == b"archived manifest"


def test_partial_download_is_not_published(download_setup):
    state = download_setup

    def interrupted(staging):
        (staging / "input.vcf.gz").write_bytes(b"partial")
        raise OSError("connection lost")

    state.download = interrupted
    with pytest.raises(OSError, match="connection lost"):
        download_synapse.main()
    assert not state.destination.exists()
    assert list(state.destination.parent.iterdir()) == []


@pytest.mark.parametrize("field,value", [
    ("id", "syn456"), ("versionNumber", 3), ("dataFileHandleId", "456"),
])
def test_download_identity_mismatch_is_not_published(download_setup, field, value):
    state = download_setup
    state.entity_updates[field] = value
    with pytest.raises(RuntimeError, match="identity differs"):
        download_synapse.main()
    assert not state.destination.exists()


@pytest.mark.parametrize("field,value,error", [
    ("fileSize", 8, "size differs"), ("md5", "0" * 32, "MD5 differs"),
])
def test_corrupt_download_is_not_published(download_setup, field, value, error):
    state = download_setup
    setattr(state.metadata, field, value)
    with pytest.raises(RuntimeError, match=error):
        download_synapse.main()
    assert not state.destination.exists()


def test_trusted_hash_mismatch_is_not_published(download_setup):
    state = download_setup
    state.argv.extend(["--expected-sha256", "0" * 64])
    with pytest.raises(RuntimeError, match="SHA-256 differs"):
        download_synapse.main()
    assert not state.destination.exists()


def test_explicit_version_mismatch_stops_before_download(download_setup):
    state = download_setup
    state.argv.extend(["--version", "1"])
    with pytest.raises(RuntimeError, match="requested version"):
        download_synapse.main()
    assert len(state.calls) == 1
    assert not state.destination.exists()


@pytest.mark.parametrize("mode", ["missing", "outside", "symlink", "directory"])
def test_unusable_download_path_is_rejected(download_setup, tmp_path, mode):
    state = download_setup
    outside = tmp_path / "existing.vcf.gz"
    outside.write_bytes(b"fixture")

    def unusable(staging):
        if mode == "missing":
            return None
        if mode == "outside":
            return outside
        if mode == "directory":
            return staging
        link = staging / "input.vcf.gz"
        link.symlink_to(outside)
        return link

    state.download = unusable
    with pytest.raises(RuntimeError, match="did not return"):
        download_synapse.main()
    assert not state.destination.exists()
    assert outside.read_bytes() == b"fixture"


def test_manifest_filename_collision_is_rejected(download_setup):
    state = download_setup
    state.name = "manifest.json"
    with pytest.raises(RuntimeError, match="conflicts with manifest"):
        download_synapse.main()
    assert not state.destination.exists()


def test_concurrent_destination_is_preserved(download_setup):
    state = download_setup
    original_download = state.download

    def competing(staging):
        state.destination.mkdir()
        (state.destination / "manifest.json").write_bytes(b"other download")
        return original_download(staging)

    state.download = competing
    with pytest.raises(FileExistsError):
        download_synapse.main()
    assert list(state.destination.iterdir()) == [state.destination / "manifest.json"]
    assert (state.destination / "manifest.json").read_bytes() == b"other download"


def test_manifest_publication_failure_removes_only_new_download(download_setup, monkeypatch):
    state = download_setup
    rename = Path.rename

    def fail_manifest(self, destination):
        if self.name == "manifest.json":
            raise OSError("publication failed")
        return rename(self, destination)

    monkeypatch.setattr(Path, "rename", fail_manifest)
    with pytest.raises(OSError, match="publication failed"):
        download_synapse.main()
    assert not state.destination.exists()


@pytest.mark.parametrize("path", ["/tmp/outside", "challenges/../../outside", "outside"])
def test_registry_path_cannot_escape_data_root(download_setup, path):
    state = download_setup
    state.registry.write_text(state.registry.read_text().replace("challenges/demography/bottleneck", path))
    with pytest.raises(ValueError, match="beneath challenges"):
        download_synapse.main()
    assert not state.calls
