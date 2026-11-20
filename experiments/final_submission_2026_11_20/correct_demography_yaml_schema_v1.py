"""Format-only repair of four rejected demographic submission YAML files."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import yaml


TASKS = {
    1: ("syn74933357", "edf4ce5516d408afd268ae6d9f3295bd475e9bb932126af92453d08d34e8d553"),
    2: ("syn74933370", "bef1635f3f75ba1bcc3c5af8d18f40899341e1175de36e082883d35736efc0bb"),
    3: ("syn74933480", "ea9254ea49330527fd23ff6d7cd41f30d80dd6a3dac6e82678c91ebe47102465"),
    4: ("syn74933479", "3ee0a388789a7969044456d82cec0bad9a49c17b04d227554f0308f12a408546"),
}
TMP = Path("/private/tmp/ghist-final-submissions-20261120")
REPO = Path(__file__).resolve().parent


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def save_new(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        assert path.read_bytes() == data, f"Refusing to replace different artifact: {path}"
    else:
        path.write_bytes(data)


def main() -> None:
    generated = datetime.now(timezone.utc).isoformat()
    for task, (entity, expected_sha) in TASKS.items():
        folder = TMP / f"task{task}"
        original = folder / "ready_for_root" / f"task{task}.final.yaml"
        raw = original.read_bytes()
        assert sha(raw) == expected_sha
        preview = folder / "official_template_preview.yaml"
        preview_raw = preview.read_bytes()
        preview_text = preview_raw.decode("utf-8")
        template = yaml.safe_load(preview_text)
        old = yaml.safe_load(raw)
        assert list(old) == list(template["parameters"])
        assert preview_text.count("parameters:\n") == 1
        metadata_text, placeholders = preview_text.split("parameters:\n")
        assert all(not line or line.startswith("  ") for line in placeholders.splitlines())
        assert raw.endswith(b"\n") and all(
            line and not line.startswith(b" ") and line.count(b": ") == 1
            for line in raw.splitlines()
        )
        nested = b"".join(b"  " + line for line in raw.splitlines(keepends=True))
        corrected = metadata_text.encode("utf-8") + b"parameters:\n" + nested
        new = yaml.safe_load(corrected)
        assert list(new) == list(template)
        assert new["parameters"] == old
        assert all(new[k] == value for k, value in template.items() if k != "parameters")
        assert original.read_bytes() == raw
        name = f"task{task}_final_schema_corrected.yaml"
        provenance = {
            "schema_version": 1,
            "created_at_utc": generated,
            "task": task,
            "scope": "Format-only correction after server INVALID: missing parameters section",
            "original": {"path": str(original), "sha256": sha(raw), "bytes": len(raw)},
            "corrected": {"filename": name, "sha256": sha(corrected), "bytes": len(corrected)},
            "official_template": {
                "entity": entity,
                "version": 1,
                "url": f"https://www.synapse.org/Synapse:{entity}",
                "preview_path": str(preview),
                "preview_sha256": sha(preview_raw),
                "preview_bytes": len(preview_raw),
                "byte_identity_claim": False,
                "evidence": "Authenticated UI preview supplied by root; display text has appended newline",
                "ui_evidence_path": str(folder / "official_template.ax.txt"),
                "ui_evidence_sha256": sha((folder / "official_template.ax.txt").read_bytes()),
            },
            "operations": [
                "Copy official descriptive/unit metadata verbatim from authenticated preview",
                "Insert original flat parameter lines under parameters with two-space indentation",
            ],
            "parameter_names_order_and_numeric_lexemes_unchanged": True,
            "original_prediction_values_unchanged": True,
            "original_preserved": True,
            "model_refit_selection_or_truth_access": False,
            "scientific_qualification_changed": False,
            "server_validation_claim": False,
            "source_sha256": sha(Path(__file__).read_bytes()),
        }
        pdata = (json.dumps(provenance, indent=2) + "\n").encode()
        for destination in (folder, REPO / f"task{task}"):
            save_new(destination / name, corrected)
            save_new(destination / f"task{task}_final_schema_correction_provenance.json", pdata)
            if destination != folder:
                save_new(destination / "official_template_preview.yaml", preview_raw)
                save_new(destination / "official_template.ax.txt", (folder / "official_template.ax.txt").read_bytes())
        print(json.dumps({"task": task, "filename": name, "sha256": sha(corrected)}))


if __name__ == "__main__":
    main()
