"""Descriptive error strata from retained raw rows; never candidate selection."""
import hashlib
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


def key(row):
    return tuple(sorted((row["indID1"], row["indID2"])))


root = Path(sys.argv[1]).resolve()
plan = load(root / "crossed_development_plan_v1.json")
report_path = root / "crossed_context_fit_v1/report.json"
report = load(report_path)
out = root / "crossed_error_strata_v1.json"
if out.exists():
    raise FileExistsError(out)
counts = defaultdict(Counter)
correct = defaultdict(Counter)
confusions = defaultdict(Counter)
bounds = {str(report_path): sha(report_path)}
for split in ("calibration", "development"):
    for role in ("base", "context"):
        scoring = root / "crossed_context_fit_v1" / f"scoring_{split}_{role}.json"
        if sha(scoring) != report["outputs_sha256"][scoring.name]:
            raise ValueError("Scoring input changed")
        bounds[str(scoring)] = sha(scoring)
        for index, cohort in enumerate(load(scoring)["cohorts"]):
            directory = root / "results_crossed_v1" / f"{split}_{index:03d}"
            result, truth = load(directory / "result.json"), load(directory / "truth.json")
            if sha(directory / "truth.json") != result["output_sha256"]["truth.json"]:
                raise ValueError("Retained truth changed")
            if truth["pairs"] != cohort["truth"]:
                raise ValueError("Scorer truth mismatch")
            group_path = root / "crossed_context_fit_v1" / f"groups_{split}_{index:02d}.json"
            if sha(group_path) != report["outputs_sha256"][group_path.name]:
                raise ValueError("Grouping audit changed")
            groups = load(group_path)["predicted_groups"]
            group_size = {sample: len(group) for group in groups for sample in group}
            config = result["cohort"]
            predictions = {key(row): row for row in cohort["predictions"]}
            if len(predictions) != len(cohort["predictions"]):
                raise ValueError("Duplicate prediction")
            for row in cohort["truth"]:
                pair = key(row)
                family = row["family"]
                ancestry = truth["pedigrees"][family]["kind"]
                hbd = max(truth["true_hbd_cm"][s] / truth["genome_cm"] for s in pair)
                hbd_bin = "under005" if hbd < .05 else "005_to015" if hbd < .15 else "atleast015"
                sizes = sorted(group_size.get(s, 1) for s in pair)
                context_bin = "both_singleton" if sizes[1] == 1 else "one_singleton" if sizes[0] == 1 else "both_grouped"
                strata = {"overall": "all", "population_size": str(config["population_size"]),
                    "missing_rate": str(config["missing_rate"]), "ancestry": ancestry,
                    "topology": "mixed" if config["complex_topology"] else "simple", "true_hbd_max": hbd_bin,
                    "predicted_sibling_context": context_bin,
                    "population_by_ancestry": f"{config['population_size']}:{ancestry}"}
                predicted = predictions.get(pair, {}).get("relation", "MISSING")
                for dimension, stratum in strata.items():
                    bucket = (split, role, dimension, stratum)
                    counts[bucket][row["relation"]] += 1
                    correct[bucket][row["relation"]] += int(predicted == row["relation"])
                    confusions[bucket][(row["relation"], predicted)] += 1
            bounds[str(directory / "result.json")] = sha(directory / "result.json")
            bounds[str(directory / "truth.json")] = sha(directory / "truth.json")
            bounds[str(group_path)] = sha(group_path)
results = []
for bucket, total in sorted(counts.items()):
    split, role, dimension, stratum = bucket
    n, n_correct = sum(total.values()), sum(correct[bucket].values())
    results.append({"split": split, "model": role, "dimension": dimension, "stratum": stratum,
        "pairs": n, "correct": n_correct, "accuracy": n_correct / n,
        "per_label": {label: {"n": total[label], "correct": correct[bucket][label],
            "recall": correct[bucket][label] / total[label] if total[label] else None} for label in plan["labels"]},
        "confusion": [{"truth": a, "prediction": b, "count": n} for (a, b), n in sorted(confusions[bucket].items())]})
for path, expected in bounds.items():
    if sha(path) != expected:
        raise ValueError("Diagnostic input changed")
output = {"scope": "Descriptive raw-row strata; no parameter/model selection, no new validation or gate changes",
    "strata": results, "bound_inputs": bounds, "script_sha256": sha(__file__), "submission_ready": False}
with out.open("x") as stream:
    stream.write(json.dumps(output, indent=2, allow_nan=False) + "\n")
print(json.dumps({"output": str(out), "sha256": sha(out), "strata": len(results)}))
