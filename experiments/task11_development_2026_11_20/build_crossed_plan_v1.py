"""Freeze a fully crossed216-cohort schedule and untouched future seeds."""
import ast
import hashlib
import itertools
import json
from collections import Counter
from pathlib import Path

root = Path(__file__).resolve().parent
remote = Path("/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task11-development-20261120")
base = json.loads((root / "protocol_v5.json").read_text())
labels = base["labels"]
kinds = ("unrelated", "first_cousins", "siblings")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def bit(text):
    return bool(hashlib.sha256(text.encode()).digest()[0] & 1)


cohorts = []
for split, replicas, seed_base in (("train", 2, 2026123000), ("calibration", 1, 2026124000),
                                   ("development", 1, 2026125000)):
    for index, (replica, ne_index, missing_index, shift, topology) in enumerate(itertools.product(range(replicas), range(3), range(3), range(3), range(2))):
        seed = seed_base + index
        third = []
        for label in labels:
            key = f"crossed-v1:train:{ne_index}:{missing_index}:{shift}:{topology}:{label}" if split == "train" else f"crossed-v1:{split}:{seed}:{label}"
            choice = bit(key) ^ bool(replica) if split == "train" else bit(key)
            third.append(True if topology and label.startswith("avuncular") else choice)
        config = {"split": split, "index": index, "global_index": len(cohorts), "seed": seed,
            "population_size": (500, 1000, 2000)[ne_index], "missing_rate": (0., .005, .02)[missing_index],
            "min_branch_size": 1, "max_branch_size": 4, "complex_topology": bool(topology),
            "ancestry_shift": shift, "replica": replica,
            "ancestry_kinds": [kinds[(label_index + shift) % 3] for label_index in range(8)],
            "third_branches": third}
        if split == "train":
            config["oof_fold"] = (2 * replica + topology + ne_index + missing_index + shift) % 4
        cohorts.append(config)
reserved = {}
for split, seed_base in (("stress", 2026126000), ("confirmation", 2026127000)):
    assignments = []
    for label in labels:
        sequence = sorted(range(8), key=lambda i: hashlib.sha256(f"crossed-v1:{split}:ancestry:{label}:{i}".encode()).digest())
        assignments.append({i: kinds[j % 3] for j, i in enumerate(sequence)})
    reserved[split] = []
    for index in range(8):
        seed = seed_base + index
        topology = index >= 4
        reserved[split].append({"split": split, "index": index, "seed": seed,
            "population_size": (500, 1000, 2000, 500, 1000, 2000, 500, 2000)[index],
            "missing_rate": (0., .005, .02, .005, .02, 0., .005, .02)[index],
            "min_branch_size": 1, "max_branch_size": 4, "complex_topology": topology,
            "ancestry_kinds": [assignment[index] for assignment in assignments],
            "third_branches": [True if topology and label.startswith("avuncular") else bit(f"crossed-v1:{split}:{seed}:{label}") for label in labels]})
registry = json.loads((root / "sparse_context_stress_design_v1.json").read_text())
old_seeds = set(registry["existing_and_reserved_seed_union"]) | {c["seed"] for c in registry["cohorts"]}
new_seeds = [c["seed"] for c in cohorts] + [c["seed"] for configs in reserved.values() for c in configs]
assert len(new_seeds) == len(set(new_seeds)) == 232 and not old_seeds & set(new_seeds)
folds = []
for fold in range(4):
    selected = [c for c in cohorts if c["split"] == "train" and c["oof_fold"] == fold]
    folds.append({"fold": fold, "cohort_count": len(selected), "cohort_indices": [c["index"] for c in selected],
        "population_size": dict(Counter(c["population_size"] for c in selected)),
        "missing_rate": dict(Counter(c["missing_rate"] for c in selected)),
        "ancestry_shift": dict(Counter(c["ancestry_shift"] for c in selected)),
        "complex_topology": dict(Counter(c["complex_topology"] for c in selected)),
        "cells": [{key: c[key] for key in ("population_size", "missing_rate", "ancestry_shift", "complex_topology", "replica")} for c in selected]})
helper = root / "diagnose_family_context.py"
tree = ast.parse(helper.read_text())
parameters = ast.literal_eval(next(n.value for n in tree.body if isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id == "PARAMETERS" for t in n.targets)))
files = ["protocol_v5.json", "source_snapshot_v5.tar.gz", "crossed_generation_v1.sbatch", "crossed_fit_v1.sbatch", "build_crossed_plan_v1.py"]
files += [str(path.relative_to(root)) for path in sorted((root / "crossed_source_v1").glob("*.py"))]
bound = {str(remote / name): sha(root / name) for name in files}
bound[str(remote / "context_diagnostic_source_v1/diagnose_family_context.py")] = sha(helper)
plan = {"schema": "task11-crossed-development-v1", "cohorts": cohorts, "reserved_future_cohorts": reserved,
    "counts": {"train": 108, "calibration": 54, "development": 54}, "labels": labels, "gates": base["gates"],
    "context_parameters": parameters, "context_metadata_override": {
        "scope": "108 fresh training,54 fresh calibration,54 fresh development; all prior cohorts excluded from fitting/calibration",
        "fold_assignment": "Explicit balanced oof_fold in each config; supersedes legacy helper metadata saying index modulo four"},
    "fold_cross_tabs": folds, "bound_inputs": bound,
    "fit_script_sha256": sha(root / "crossed_source_v1/fit_crossed_context.py"),
    "source_change": "Only genotype generator nuisance assignment changes: explicit ancestry/third-branch lists replace cohort-index coupling. Founder, gene-drop, detector,149 pair features,36 context features and classifier unchanged.",
    "design": "Train two complete3Ne×3missing×3ancestry-shift×2topology blocks; calibration/development one each. Branch sizes1–4 independently sampled by each master genotype seed. All216 outcomes retained.",
    "metadata_randomness": "SHA256 domain-separated branch draws, complemented between train replicas except required mixed avuncular branches. Future ancestry lists are independently shuffled per label and balanced3/3/2 over eight cohorts.",
    "launch_gates": "All nine original gates on54 development cohorts plus calibration min-label recall>=.60/ECE<=.10. Then all nine gates on fresh stress, then all nine on fresh confirmation. No parameter changes, refits or exclusions after development.",
    "bootstrap": {"replicates": 2000, "seed": 20301120, "unit": "whole cohort"},
    "old_seed_union": sorted(old_seeds), "new_seed_union": sorted(new_seeds), "no_seed_overlap": True,
    "old_registry_sha256": sha(root / "sparse_context_stress_design_v1.json"),
    "budget": {"partition": "genoa-std-mem", "max_concurrent_generation": 4, "cpus_per_task": 4,
        "mem_per_cpu_mb": 5632, "generation_walltime_minutes": 10, "fit_walltime_minutes": 25,
        "latest_generation_start_utc": "2026-11-20T21:40:00+00:00", "latest_fit_start_utc": "2026-11-20T21:50:00+00:00",
        "development_decision_due_utc": "2026-11-20T22:15:00+00:00", "latest_stress_confirmation_start_utc": "2026-11-20T22:25:00+00:00",
        "handoff_due_utc": "2026-11-20T22:45:00+00:00", "estimated_generation_minutes": [80, 120]},
    "preservation_review": "Report base-correct flips and high-support>=.9 flips; no outcome-selected guard change.",
    "limitations": ["Finite synthetic pedigree templates and coalescent founders remain transfer assumptions.",
                    "Earlier failed candidates and their evidence remain immutable; no previous development cohort is corrected confirmation."],
    "submission_ready": False, "final_challenge_access": False, "submission_count": 0}
path = root / "crossed_development_plan_v1.json"
with path.open("x") as destination:
    json.dump(plan, destination, indent=2)
    destination.write("\n")
print(sha(path))
