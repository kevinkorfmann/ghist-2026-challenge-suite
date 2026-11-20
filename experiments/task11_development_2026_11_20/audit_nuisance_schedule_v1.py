"""Verify actual retained ancestry kinds against population-size schedules."""
import hashlib
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

root = Path(sys.argv[1]).resolve()
protocol = json.loads((root / "protocol_v5.json").read_text())
labels, kinds = protocol["labels"], ("unrelated", "first_cousins", "siblings")
cells, inputs = defaultdict(Counter), {}
rows = []
for split, indices in (("old_train", range(16)), ("new_train", range(16, 64)),
                       ("calibration", range(4)), ("exposed", range(8))):
    for index in indices:
        if split.endswith("train"):
            directory = root / "expanded_training_raw_metadata_v1" / f"train_{index:02d}"
            paths = {name: directory / f"{name}.json" for name in ("truth", "result")}
        else:
            paths = {name: root / "expanded_development_audit_v1" / f"{split}_{index:02d}_{name}.json" for name in ("truth", "result")}
        truth, result = (json.loads(paths[name].read_text()) for name in ("truth", "result"))
        config = result["cohort"]
        for path in paths.values():
            inputs[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
        for label_index, (label, pedigree) in enumerate(zip(labels, truth["pedigrees"], strict=True)):
            expected = kinds[(config["index"] + label_index) % 3]
            if pedigree["kind"] != expected:
                raise ValueError("Retained ancestry does not match source schedule")
            row = {"split": split, "cohort_seed": config["seed"], "cohort_index": config["index"],
                   "focal_label": label, "population_size": config["population_size"], "ancestry_kind": pedigree["kind"]}
            rows.append(row)
            cells[(split, label, config["population_size"])][pedigree["kind"]] += 1
cross_tab = [{"split": split, "focal_label": label, "population_size": ne, "ancestry_counts": dict(counts)}
             for (split, label, ne), counts in sorted(cells.items())]
all_deterministic = all(len(counts) == 1 for counts in cells.values())
rotated = sum(set(cells[("old_train", label, ne)]) != set(cells[("new_train", label, ne)])
              for label in labels for ne in (500, 1000, 2000))
output = {"scope": "Verified design confounding; causal effect on classification remains a hypothesis",
    "retained_cohorts": len(rows) // 8, "retained_families": len(rows), "all_recorded_kinds_match_source": True,
    "each_split_label_population_cell_has_one_ancestry_kind": all_deterministic,
    "old_to_new_label_population_cells_rotated": rotated, "total_label_population_cells": 24,
    "cross_tab": cross_tab, "rows": rows, "inputs_sha256": inputs,
    "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "interpretation": "Ne and ancestry kind are perfectly coupled conditional on focal label in each old/new block. All24 label×Ne cells rotate in added training. Corrected training must cross or independently randomize nuisance assignments. This audit does not prove how much regression this caused."}
path = root / "nuisance_schedule_audit_v1.json"
with path.open("x") as destination:
    json.dump(output, destination, indent=2)
    destination.write("\n")
print(json.dumps({key: value for key, value in output.items() if key not in ("rows", "inputs_sha256", "cross_tab")}, indent=2))
