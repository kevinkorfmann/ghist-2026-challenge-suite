import json
from types import SimpleNamespace

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
                return SimpleNamespace(versionNumber=2)
            destination = tmp_path / kwargs["downloadLocation"] / "input.vcf.gz"
            destination.write_bytes(b"fixture")
            return SimpleNamespace(path=str(destination), dataFileHandleId="123")

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
