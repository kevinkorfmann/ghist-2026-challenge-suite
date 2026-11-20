"""Independent official-template and exact-scalar audit; imports no repair code."""

from __future__ import annotations

import hashlib
import json
import re
from datetime import datetime, timezone
from decimal import Decimal
from pathlib import Path

import yaml


ROOT = Path("/private/tmp/ghist-final-submissions-20261120")
REPO = Path(__file__).resolve().parent
KEYS = {
    1: ["generations", "post_decline_fraction"],
    2: ["mainland_population_size", "island_population_size", "generations_since_split",
        "generations_since_migration", "migration_rate"],
    3: ["generations", "post_growth_scaling"],
    4: ["AFR_population_size", "EUR_population_size", "generations_since_split", "migration_rate"],
}
ENTITIES = {1: "syn74933357", 2: "syn74933370", 3: "syn74933480", 4: "syn74933479"}


def digest(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def scalar_map(node: yaml.Node) -> dict[str, tuple[str, str]]:
    assert isinstance(node, yaml.MappingNode)
    result = {}
    for key, value in node.value:
        assert isinstance(key, yaml.ScalarNode) and isinstance(value, yaml.ScalarNode)
        assert key.value not in result
        result[key.value] = (value.tag, value.value)
    return result


def main() -> None:
    for task, keys in KEYS.items():
        folder = ROOT / f"task{task}"
        original_path = folder / "ready_for_root" / f"task{task}.final.yaml"
        corrected_path = folder / f"task{task}_final_schema_corrected.yaml"
        original = original_path.read_bytes()
        corrected = corrected_path.read_bytes()
        preview = (folder / "official_template_preview.yaml").read_bytes()
        ax = (folder / "official_template.ax.txt").read_text()
        assert ENTITIES[task] in ax.splitlines()[0]
        matches = re.findall(r"Size (\d+) Bytes MD5 ([0-9a-f]{32})", ax)
        assert len(matches) == 1
        official_size, official_md5 = int(matches[0][0]), matches[0][1]
        template_candidates = [preview]
        if preview.endswith(b"\n"):
            template_candidates.append(preview[:-1])
        exact = [candidate for candidate in template_candidates
                 if len(candidate) == official_size
                 and hashlib.md5(candidate).hexdigest() == official_md5]
        assert len(exact) == 1, "Preview must recover exact official bytes after optional added newline"
        official_bytes = exact[0]
        template = yaml.safe_load(official_bytes)
        assert list(template["parameters"]) == keys
        old_node = yaml.compose(original)
        old_scalars = scalar_map(old_node)
        assert list(old_scalars) == keys
        new_node = yaml.compose(corrected)
        assert isinstance(new_node, yaml.MappingNode)
        top_keys = [key.value for key, _ in new_node.value]
        assert len(top_keys) == len(set(top_keys)) and top_keys == list(template)
        parameters = [value for key, value in new_node.value if key.value == "parameters"]
        assert len(parameters) == 1
        new_scalars = scalar_map(parameters[0])
        assert old_scalars == new_scalars and list(new_scalars) == keys
        for tag, value in new_scalars.values():
            assert tag in {"tag:yaml.org,2002:int", "tag:yaml.org,2002:float"}
            assert Decimal(value).is_finite() and Decimal(value) > 0
        new = yaml.safe_load(corrected)
        assert new["parameters"] == yaml.safe_load(original)
        assert all(new[key] == value for key, value in template.items() if key != "parameters")
        # Independent line comparison proves no float serializer changed any digit.
        parameter_lines = corrected.split(b"parameters:\n", 1)[1].splitlines(keepends=True)
        assert b"".join(line[2:] for line in parameter_lines) == original
        assert original_path.read_bytes() == original
        assert (REPO / f"task{task}" / corrected_path.name).read_bytes() == corrected
        provenance_path = folder / f"task{task}_final_schema_correction_provenance.json"
        provenance = json.loads(provenance_path.read_text())
        assert provenance["original"]["sha256"] == digest(original)
        assert provenance["corrected"]["sha256"] == digest(corrected)
        assert (REPO / f"task{task}" / provenance_path.name).read_bytes() == provenance_path.read_bytes()
        report = {
            "task": task, "status": "PASS",
            "checked_at_utc": datetime.now(timezone.utc).isoformat(),
            "original_sha256": digest(original), "corrected_sha256": digest(corrected),
            "official_template_entity": ENTITIES[task], "official_template_version": 1,
            "official_template_md5": official_md5, "official_template_bytes": official_size,
            "official_template_exact_bytes_sha256": digest(official_bytes),
            "exact_template_recovered_from_preview": True,
            "preview_appended_newline_removed": official_bytes != preview,
            "independent_hardcoded_task_keys_match_template": True,
            "parameter_names_order_numeric_tags_and_lexemes_exact": True,
            "complete_original_lines_recovered_after_unindent": True,
            "metadata_matches_template": True, "finite_positive_parameters": True,
            "duplicate_keys": False, "original_unchanged": True,
            "temporary_and_repository_copies_identical": True,
            "model_refit_or_selection_changes": False,
            "scientific_qualification_changed": False, "server_validated": False,
            "verifier_sha256": digest(Path(__file__).read_bytes()),
            "verifier_imports_repair_code": False,
        }
        payload = (json.dumps(report, indent=2) + "\n").encode()
        for destination in (folder, REPO / f"task{task}"):
            path = destination / f"task{task}_final_schema_independent_check.json"
            assert not path.exists()
            path.write_bytes(payload)
        print(json.dumps({"task": task, "status": "PASS", "sha256": digest(corrected),
                          "official_template_exact_sha256": digest(official_bytes)}))


if __name__ == "__main__":
    main()
