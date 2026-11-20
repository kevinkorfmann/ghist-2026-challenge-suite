"""Recompute fixed development launch gates and stage compact raw evidence."""
import datetime
import hashlib
import importlib.util
import json
import shutil
import sys
from collections import Counter
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


root = Path(sys.argv[1]).resolve()
helper = root / "context_diagnostic_source_v1/diagnose_family_context.py"
if sha(helper) != "4b347c189e0a5799c99e158715c236934485e618db44cb231e84a27535fd396b":
    raise ValueError("Context source changed")
spec = importlib.util.spec_from_file_location("frozen_context", helper)
ctx = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ctx)
dev, np = ctx.dev, ctx.np
protocol = dev.load(root / "protocol_v5.json")
gates = dev.validate_protocol(protocol)
dev.require_hashes({str(dev.ROOT / path): value for path, value in protocol["source_sha256"].items()})
fit = root / "crossed_context_fit_v1"
report = dev.load(fit / "report.json")
dev.require_hashes(report["inputs_sha256"])
dev.require_hashes({str(fit / path): value for path, value in report["outputs_sha256"].items()})
design_path = root / "crossed_development_plan_v1.json"
dev.require_hashes({str(design_path): "e90d2235de805d17a6c1d0176cab2132c2626b5ef04d3c97b4a4e12c67850016"})
plan = dev.load(design_path)
if plan["gates"] != gates:
    raise ValueError("Crossed gates changed")
out = root / "crossed_development_audit_v1"
out.mkdir(exist_ok=False)
summaries, partitions = {}, {}
for split, count in (("calibration", 54), ("development", 54)):
    partitions[split] = []
    for index in range(count):
        directory = root / "results_crossed_v1" / f"{split}_{index:03d}"
        result = dev.load(directory / "result.json")
        for filename in ("truth.json", "clusters.json", "protocol.json"):
            dev.require_hashes({str(directory / filename): result["output_sha256"][filename]})
            shutil.copyfile(directory / filename, out / f"{split}_{index:02d}_{filename}")
        shutil.copyfile(directory / "result.json", out / f"{split}_{index:02d}_result.json")
        truth = dev.load(directory / "truth.json")
        partition = dev.partition_metrics(truth["families"], dev.load(directory / "clusters.json")["families"])
        if partition != result["partition"]:
            raise ValueError("Stored partition differs from raw families")
        partitions[split].append(partition)
    for model in ("base", "context"):
        filename = f"scoring_{split}_{model}.json"
        scorer = dev.load(fit / filename)
        if scorer["labels"] != list(dev.LABELS) or len(scorer["cohorts"]) != count:
            raise ValueError("Incomplete raw scoring input")
        truth, rows, scores = {"pairs": []}, [], []
        for index, cohort in enumerate(scorer["cohorts"]):
            if cohort["cohort_id"] != f"{split}_{index:02d}":
                raise ValueError("Cohort order mismatch")
            retained = dev.load(out / f"{split}_{index:02d}_truth.json")
            if retained["pairs"] != cohort["truth"]:
                raise ValueError("Raw scorer truth differs from generator truth")
            for row in cohort["predictions"]:
                probabilities = row["probabilities"]
                if set(probabilities) != set(dev.LABELS) or abs(sum(probabilities.values()) - 1) > 1e-12:
                    raise ValueError("Invalid probabilities")
                if any(not np.isfinite(p) or p < 0 for p in probabilities.values()):
                    raise ValueError("Invalid probability value")
                if abs(max(probabilities.values()) - row["support"]) > 1e-12:
                    raise ValueError("Stored support mismatch")
                if probabilities[row["relation"]] != max(probabilities.values()):
                    raise ValueError("Stored prediction differs from probabilities")
            scores.append(dev.score({"pairs": cohort["truth"]}, cohort["predictions"]))
            for source, target in ((cohort["truth"], truth["pairs"]), (cohort["predictions"], rows)):
                target.extend({**row, "indID1": f"c{index}_{row['indID1']}", "indID2": f"c{index}_{row['indID2']}"} for row in source)
        overall = dev.score(truth, rows)
        if overall != report["reports"][split][model]:
            raise ValueError("Stored fit summary differs from raw truth/predictions")
        rng = np.random.default_rng(20301120)
        bootstrap = []
        for _ in range(2000):
            correct, counts = Counter(), Counter()
            for index in rng.integers(0, count, count):
                correct.update(scores[index]["correct_counts"])
                counts.update(scores[index]["counts"])
            bootstrap.append(np.mean([correct[label] / counts[label] for label in dev.LABELS]))
        overall["cohort_bootstrap_balanced_interval95"] = np.quantile(bootstrap, [.025, .975]).tolist()
        exact = sum(p["exactly_recovered_families"] for p in partitions[split]) / (8 * count)
        missing = sum(len(p["missing_samples"]) for p in partitions[split])
        checks = {"accuracy": overall["accuracy"] >= gates["accuracy_min"],
            "balanced_accuracy": overall["balanced_accuracy"] >= gates["balanced_accuracy_min"],
            "each_label_recall": min(overall["recall"].values()) >= gates["each_label_recall_min"],
            "bootstrap_lower": overall["cohort_bootstrap_balanced_interval95"][0] >= gates["bootstrap_balanced_lower95_min"],
            "family_recovery": exact >= gates["exact_family_recovery_min"],
            "missing_samples": missing <= gates["missing_samples_max"],
            "missing_true_pairs": overall["missing_true_pairs"] <= gates["missing_true_pairs_max"],
            "extra_pairs": overall["extra_cross_family_pairs"] <= gates["extra_cross_family_pairs_max"],
            "calibration": overall["ece_10_equal_width"] <= gates["expected_calibration_error_max"]}
        summaries[f"{split}_{model}"] = {"score": overall, "checks": checks,
            "exact_family_recovery": exact, "all_nine_gates_passed": all(checks.values())}
        shutil.copyfile(fit / filename, out / filename)
calibration = summaries["calibration_context"]["checks"]
scientific_launch = summaries["development_context"]["all_nine_gates_passed"] and calibration["each_label_recall"] and calibration["calibration"]
within_budget = datetime.datetime.now(datetime.UTC) <= datetime.datetime.fromisoformat(plan["budget"]["development_decision_due_utc"])
launch = scientific_launch and within_budget
for path in (design_path, fit / "report.json", root / "protocol_v5.json"):
    shutil.copyfile(path, out / path.name)
summary = {"scope": "Development launch guard, not scientific qualification", "summaries": summaries,
    "fresh_stress_launch_allowed": launch, "scientific_development_gates_passed": scientific_launch,
    "within_decision_budget": within_budget, "calibration_gate_scope": ["each_label_recall", "calibration"],
    "gates": gates, "preservation_audit": {split: report["reports"][split]["flips"] for split in ("calibration", "development")},
    "models_sha256": {name: sha(fit / name) for name in ("base_model.joblib", "context_model.joblib")},
    "fit_report_sha256": sha(fit / "report.json"), "audit_script_sha256": sha(__file__),
    "input_design_sha256": sha(design_path), "launch_gate_policy": plan["launch_gates"],
    "output_sha256": {path.name: sha(path) for path in sorted(out.iterdir())},
    "submission_ready": False, "reserved_confirmation_used": False}
dev.dump(out / "audit.json", summary)
print(json.dumps({"summaries": summaries, "fresh_stress_launch_allowed": launch}, indent=2))
