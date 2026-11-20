"""One prospectively fixed 64-cohort development fit and matched context ablation."""
import hashlib
import importlib.util
import json
import os
import sys
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


root = Path(sys.argv[1]).resolve()
plan_path = root / "expanded_training_plan_v2.json"
plan = json.loads(plan_path.read_text())
for path, expected in plan["bound_inputs"].items():
    if sha(path) != expected:
        raise ValueError(f"Expanded plan input changed: {path}")
if sha(__file__) != plan["fit_script_sha256"]:
    raise ValueError("Expanded fit source changed")
if "SLURM_JOB_ID" not in os.environ:
    raise RuntimeError("Compute requires Slurm")
helper = root / "context_diagnostic_source_v1/diagnose_family_context.py"
spec = importlib.util.spec_from_file_location("frozen_context", helper)
ctx = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ctx)
dev, np = ctx.dev, ctx.np
protocol = dev.load(root / "protocol_v5.json")
dev.require_runtime(protocol)
dev.require_hashes(protocol["inputs"])
dev.require_hashes({str(dev.ROOT / p): h for p, h in protocol["source_sha256"].items()})
assert ctx.PARAMETERS == plan["context_parameters"]
out = root / "expanded_context_fit_v2"
out.mkdir(exist_ok=False)
dev.dump(out / "prospective_plan.json", plan)
dev.dump(out / "job.json", {"slurm_job_id": os.environ["SLURM_JOB_ID"]})
bounds = {**plan["bound_inputs"], str(plan_path): sha(plan_path), str(Path(__file__).resolve()): sha(__file__)}
data, paths = {}, {}
generation_protocol = dev.load(root / "expanded_training_protocol_v1.json")
for split, count in (("train", 64), ("calibration", 4), ("exposed", 8)):
    data[split], paths[split] = [], []
    for index in range(count):
        if split == "train" and index >= 16:
            directory = root / "results_expanded_v1" / f"train_{index:02d}"
        else:
            directory = root / ("results_v3" if split == "exposed" else "results_v5") / f"{'confirmation' if split == 'exposed' else split}_{index:02d}"
        d = list(dev.dataset(directory / "result.json"))
        if split == "train" and index >= 16:
            assert d[0]["cohort"] == generation_protocol["cohorts"][index - 16]
            assert dev.load(directory / "protocol.json")["master_sha256"] == sha(root / "expanded_training_protocol_v1.json")
        if split == "exposed":
            phase_path = root / "phase_ablation_v5" / f"exposed_cohort_{index:02d}.npz"
            expected = dev.load(root / "phase_development_ablation_v5.json")["cohorts"][index]
            dev.require_hashes({str(phase_path): expected["phase_features_sha256"]})
            d[3] = np.concatenate([d[3], np.load(phase_path)["features"]], axis=1)
            bounds[str(phase_path)] = sha(phase_path)
        data[split].append(d)
        paths[split].append(directory)
        bounds[str(directory / "result.json")] = sha(directory / "result.json")
assert len({d[0]["cohort"]["seed"] for cohorts in data.values() for d in cohorts}) == 76
assert sum(d[0]["cohort"]["complex_topology"] for d in data["train"]) == 32


def temperature(model, calibration):
    x = np.concatenate([d[3][d[4]] for d in calibration])
    y = np.concatenate([d[5] for d in calibration])
    losses = {str(t): float(dev.log_loss(y, dev.temperature_scale(model.predict_proba(x), t), labels=model.classes_))
              for t in ctx.PARAMETERS["temperature_grid"]}
    return float(min(losses, key=losses.get)), losses


base_model = ctx.train_model(data["train"])
dev.joblib.dump(base_model, out / "base_model.joblib")
base_temperature, base_losses = temperature(base_model, data["calibration"])
fold_models = {}
for fold in range(4):
    fold_models[fold] = ctx.train_model([d for i, d in enumerate(data["train"]) if i % 4 != fold])
    dev.joblib.dump(fold_models[fold], out / f"sibling_fold_{fold}.joblib")
map_path = root.parents[1] / "external/dog-recombination/campbell2016_average.plink.map"
maps = dev.read_map(map_path)
augmented, base_predictions = {}, {}
for split, cohorts in data.items():
    augmented[split], base_predictions[split] = [], []
    for index, (d, directory) in enumerate(zip(cohorts, paths[split], strict=True)):
        base = dev.predict(d[2], d[3], base_model, base_temperature, False)
        group_model = fold_models[index % 4] if split == "train" else base_model
        grouping = dev.predict(d[2], d[3], group_model, base_temperature, False)
        samples = {sample for row in d[2] for sample in dev.key(row)}
        records, _ = dev.read_haplotype_segments(directory / "segments.ibd.gz", samples)
        masks = {t: dev.read_hbd_masks(directory / "segments.hbd.gz", t) for t in ctx.PARAMETERS["hbd_cutoffs_cm"]}
        features, audit = ctx.context_features(d[2], ctx.probabilities(grouping), list(dev.LABELS).index("siblings"), records, masks, maps)
        membership = {sample: frozenset(group) for group in audit["predicted_groups"] for sample in group}
        actual = {dev.key(row): row["relation"] for row in d[1]["pairs"]}
        audit["grouping_truth_audit"] = {"false_grouped_pairs": sum(membership.get(a, frozenset((a,))) == membership.get(b, frozenset((b,))) and label != "siblings" for (a, b), label in actual.items()),
            "missed_sibling_pairs": sum(membership.get(a, frozenset((a,))) != membership.get(b, frozenset((b,))) and label == "siblings" for (a, b), label in actual.items()),
            "out_of_fold_grouping": split == "train"}
        np.savez_compressed(out / f"context_{split}_{index:02d}.npz", features=features)
        dev.dump(out / f"groups_{split}_{index:02d}.json", audit)
        augmented[split].append([*d[:3], np.concatenate([d[3], features], axis=1), *d[4:]])
        base_predictions[split].append(base)
context_model = ctx.train_model(augmented["train"])
dev.joblib.dump(context_model, out / "context_model.joblib")
context_temperature, context_losses = temperature(context_model, augmented["calibration"])
reports = {}
for split in ("calibration", "exposed"):
    cohorts = augmented[split]
    context_rows = [dev.predict(d[2], d[3], context_model, context_temperature, False) for d in cohorts]
    flips = {"base_correct_to_wrong": 0, "base_wrong_to_correct": 0, "strong_base_correct_to_wrong": 0}
    for index, (d, base, context) in enumerate(zip(cohorts, base_predictions[split], context_rows, strict=True)):
        actual = {dev.key(row): row["relation"] for row in d[1]["pairs"]}
        for b, c in zip(base, context, strict=True):
            label = actual.get(dev.key(b))
            flips["base_correct_to_wrong"] += int(b["relation"] == label and c["relation"] != label)
            flips["base_wrong_to_correct"] += int(b["relation"] != label and c["relation"] == label)
            flips["strong_base_correct_to_wrong"] += int(b["support"] >= .9 and b["relation"] == label and c["relation"] != label)
    reports[split] = {"base": ctx.aggregate(cohorts, base_predictions[split]),
                      "context": ctx.aggregate(cohorts, context_rows), "flips": flips}
    for name, predicted in (("base", base_predictions[split]), ("context", context_rows)):
        dev.dump(out / f"scoring_{split}_{name}.json", {"labels": list(dev.LABELS), "cohorts": [
            {"cohort_id": f"{split}_{i:02d}", "truth": d[1]["pairs"], "predictions": rows}
            for i, (d, rows) in enumerate(zip(cohorts, predicted, strict=True))]})
dev.require_hashes(bounds)
training_audit = {"cohorts": len(data["train"]), "families": 8 * len(data["train"]),
    "observed_feature_rows": sum(len(d[2]) for d in data["train"]), "fitted_true_pair_rows": sum(len(d[4]) for d in data["train"]),
    "exact_families": sum(d[0]["partition"]["exactly_recovered_families"] for d in data["train"]),
    "independent_focal_families_per_label": 64, "simple_cohorts": 32, "mixed_generation_cohorts": 32}
report = {"scope": "Expanded training development experiment; not scientific qualification or final inference",
    "training_audit": training_audit, "base_temperature": base_temperature, "context_temperature": context_temperature,
    "base_calibration_losses": base_losses, "context_calibration_losses": context_losses,
    "reports": reports, "inputs_sha256": bounds, "outputs_sha256": {p.name: sha(p) for p in sorted(out.iterdir())},
    "all_original_gates_unchanged": protocol["gates"], "submission_ready": False,
    "limitations": ["Earlier calibration/exposed cohorts have informed development and are not new confirmation.",
                    "More independent families still come from finite sibling-rich templates and coalescent population assumptions.",
                    "Context model is a matched diagnostic; preservation of strong predictions requires separate review."]}
dev.dump(out / "report.json", report)
print(json.dumps({"training": training_audit, "reports": reports}, indent=2))
