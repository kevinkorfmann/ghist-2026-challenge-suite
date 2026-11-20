"""Fixed development-only sibling-union/IBD011 diagnostic; never qualification.

Features use observed partitions, predicted sibling groups and detector intervals.
Training sibling groups come from four cohort-level out-of-fold pair models.
No truth grouping, sample identifier or pedigree property enters a feature.
"""
import json
import sys
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import numpy as np

ROOT = Path(sys.argv[1]) if len(sys.argv) > 1 else None
if ROOT is not None and (ROOT / "source_v5").is_dir():
    SOURCE = ROOT / "source_v5"
    sys.path[:0] = [str(SOURCE), str(SOURCE / "src"), str(SOURCE / "scripts")]
    from scripts import develop_task11_pedigree as dev
else:
    dev = None

PARAMETERS = {"version": "task11-context-information-diagnostic-v1", "sibling_probability_min": .95,
    "cohort_folds": 4, "fold_assignment": "training cohort index modulo four", "phase_min_cm": 10.,
    "hbd_cutoffs_cm": [2., 5., 10.], "ambiguous_component_policy": "all members become singletons",
    "edge_filter": "merge each homolog edge then require at least10cM; union eligible edges",
    "ibd0_exclusion": "all detected IBD between siblings, including short edges",
    "mask": "union of all participating individuals' HBD at each fixed cutoff",
    "zero_exposure": "IBD011 ratio, numerator and denominator are zero; union remains measurable",
    "model": {"n_estimators": 800, "min_samples_leaf": 3, "max_features": .75, "n_jobs": 4, "random_state": 20300111},
    "temperature_grid": [.75, 1., 1.25, 1.5, 2.],
    "selection": "single fixed diagnostic; all results reported, no outcome-selected feature/model variant",
    "scope": "training only fitting, four separate calibration cohorts, exposed8 development only",
    "submission_ready": False}


def merge(intervals):
    result = []
    for start, end in sorted(intervals):
        if end <= start:
            raise ValueError("Nonpositive interval")
        if result and start <= result[-1][1]:
            result[-1] = (result[-1][0], max(end, result[-1][1]))
        else:
            result.append((start, end))
    return result


def intersect(a, b):
    a, b = merge(a), merge(b)
    result, i, j = [], 0, 0
    while i < len(a) and j < len(b):
        start, end = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if start < end:
            result.append((start, end))
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return result


def subtract(intervals, excluded):
    result = []
    for start, end in merge(intervals):
        position = start
        for left, right in merge(excluded):
            if right <= position:
                continue
            if left >= end:
                break
            if left > position:
                result.append((position, min(left, end)))
            position = max(position, right)
        if position < end:
            result.append((position, end))
    return result


def length(intervals, chromosome, maps):
    bp, cm = maps[chromosome]
    total = 0.
    for start, end in merge(intervals):
        if start < bp[0] or end > bp[-1]:
            raise ValueError("Context interval exceeds map")
        total += float(np.diff(np.interp([start, end], bp, cm))[0])
    return total


def predicted_groups(rows, probabilities, sibling_index):
    """Connected threshold components accepted only if cliques, independent of IDs."""
    families, lookup = defaultdict(set), {}
    for row, probability in zip(rows, probabilities, strict=True):
        a, b = row["indID1"], row["indID2"]
        families[row["family"]].update((a, b))
        lookup[frozenset((a, b))] = float(probability[sibling_index])
    groups, audit = {}, {"ambiguous_components": 0, "predicted_groups": []}
    for members in families.values():
        if any(frozenset(pair) not in lookup for pair in combinations(members, 2)):
            raise ValueError("Incomplete observed family pairs")
        adjacency = {a: {b for b in members if a != b and lookup[frozenset((a, b))] >= .95} for a in members}
        unseen = set(members)
        while unseen:
            component, pending = set(), {next(iter(unseen))}
            while pending:
                item = pending.pop()
                if item in component:
                    continue
                component.add(item)
                pending.update(adjacency[item] - component)
            unseen -= component
            clique = all(b in adjacency[a] for a, b in combinations(component, 2))
            if not clique:
                audit["ambiguous_components"] += 1
            accepted = [component] if clique else [{a} for a in component]
            for group in accepted:
                audit["predicted_groups"].append(sorted(group))
                groups.update({a: frozenset(group) for a in group})
    audit["predicted_groups"].sort()
    return groups, audit


def interval_index(records, maps):
    all_ibd, eligible = {}, {}
    for (a, b, chromosome), values in records.items():
        edge = defaultdict(list)
        for start, end, ha, hb in values:
            edge[(ha, hb)].append((start, end))
        all_ibd[(frozenset((a, b)), chromosome)] = merge([(s, e) for s, e, *_ in values])
        eligible[(frozenset((a, b)), chromosome)] = merge([interval for values in edge.values()
            for interval in merge(values) if length([interval], chromosome, maps) >= 10.])
    return all_ibd, eligible


def side_context(relative, siblings, all_ibd, eligible, masks, maps):
    siblings = frozenset(siblings) - {relative}
    if not siblings:
        raise ValueError("No opposite endpoint in predicted group")
    genome = sum(float(cm[-1] - cm[0]) for _, cm in maps.values())
    union_cm = sum_cm = union_011_cm = sum_011_cm = exposure_cm = 0.
    triples = list(combinations(sorted(siblings), 2))
    for chromosome, (bp, _) in maps.items():
        pair_spans = {}
        for sibling in siblings:
            excluded = merge(masks.get((relative, chromosome), []) + masks.get((sibling, chromosome), []))
            pair_spans[sibling] = subtract(eligible.get((frozenset((relative, sibling)), chromosome), []), excluded)
        union_cm += length([span for spans in pair_spans.values() for span in spans], chromosome, maps)
        sum_cm += sum(length(spans, chromosome, maps) for spans in pair_spans.values())
        shared_011 = []
        for c, d in triples:
            excluded = merge(all_ibd.get((frozenset((c, d)), chromosome), []) +
                             masks.get((relative, chromosome), []) + masks.get((c, chromosome), []) +
                             masks.get((d, chromosome), []))
            exposure = subtract([(int(bp[0]), int(bp[-1]))], excluded)
            spans = intersect(intersect(pair_spans[c], pair_spans[d]), exposure)
            exposure_cm += length(exposure, chromosome, maps)
            sum_011_cm += length(spans, chromosome, maps)
            shared_011.extend(spans)
        union_011_cm += length(shared_011, chromosome, maps)
    return np.array([union_cm / genome, sum_cm / len(siblings) / genome,
        union_011_cm / genome, sum_011_cm / exposure_cm if exposure_cm else 0.,
        exposure_cm / len(triples) / genome if triples else 0.,
        sum_011_cm / len(triples) / genome if triples else 0.])


def context_features(rows, probabilities, sibling_index, records, masks_by_cutoff, maps):
    groups, audit = predicted_groups(rows, probabilities, sibling_index)
    all_ibd, eligible = interval_index(records, maps)
    output, zeros = [], []
    for cutoff in PARAMETERS["hbd_cutoffs_cm"]:
        cache, values = {}, []
        for row in rows:
            a, b = row["indID1"], row["indID2"]
            sides = []
            for relative, siblings in ((a, groups[b]), (b, groups[a])):
                key = (relative, siblings)
                if key not in cache:
                    cache[key] = side_context(relative, siblings, all_ibd, eligible, masks_by_cutoff[cutoff], maps)
                sides.append(cache[key])
            # Separate low/high values per feature preserve endpoint swap invariance.
            values.append(np.sort(np.stack(sides), axis=0).T.reshape(-1))
        matrix = np.array(values)
        output.append(matrix)
        zeros.append(int((matrix[:, 8:10] == 0).sum()))
    matrix = np.concatenate(output, axis=1)
    if not np.isfinite(matrix).all():
        raise ValueError("Nonfinite context features")
    audit["zero_ibd0_exposure_endpoints"] = zeros
    audit["feature_names"] = [f"hbd{int(t)}_{feature}_{side}" for t in PARAMETERS["hbd_cutoffs_cm"]
        for feature in ("union_ibd_fraction", "mean_ibd_fraction", "union_011_fraction", "ratio_011_to_ibd0",
                        "mean_ibd0_exposure_fraction", "mean_011_fraction") for side in ("low", "high")]
    return matrix, audit


def train_model(data):
    x = np.concatenate([d[3][d[4]] for d in data])
    y = np.concatenate([d[5] for d in data])
    model = dev.ExtraTreesClassifier(**PARAMETERS["model"])
    model.fit(x, y, sample_weight=dev.compute_sample_weight("balanced", y))
    return model


def probabilities(rows):
    return np.array([[r["probabilities"][label] for label in dev.LABELS] for r in rows])


def aggregate(data, predicted):
    truth, rows = {"pairs": []}, []
    for i, (d, p) in enumerate(zip(data, predicted, strict=True)):
        for source, target in ((d[1]["pairs"], truth["pairs"]), (p, rows)):
            target.extend({**r, "indID1": f"c{i}_{r['indID1']}", "indID2": f"c{i}_{r['indID2']}"} for r in source)
    return dev.score(truth, rows)


def main(root):
    if "SLURM_JOB_ID" not in dev.os.environ:
        raise RuntimeError("Compute requires Slurm")
    plan_path = root / "context_diagnostic_plan_v1.json"
    plan = dev.load(plan_path)
    assert plan["parameters"] == PARAMETERS
    assert plan["diagnostic_script_sha256"] == dev.digest(Path(__file__))
    dev.require_hashes(plan["bound_inputs"])
    out = root / "context_diagnostic_v1"
    out.mkdir(exist_ok=False)
    dev.dump(out / "parameters_before_results.json", PARAMETERS)
    dev.dump(out / "prospective_plan.json", plan)
    protocol = dev.load(root / "protocol_v5.json")
    dev.require_runtime(protocol)
    dev.require_hashes(protocol["inputs"])
    dev.require_hashes({str(SOURCE / p): h for p, h in protocol["source_sha256"].items()})
    candidate = dev.load(root / "candidate_v5/candidate.json")
    dev.require_hashes(candidate["bound_inputs"])
    baseline_model = dev.joblib.load(candidate["model_path"])
    split_data, split_paths = {}, {}
    bounds = {str(root / "protocol_v5.json"): dev.digest(root / "protocol_v5.json"),
              str(root / "candidate_v5/candidate.json"): dev.digest(root / "candidate_v5/candidate.json"),
              str(plan_path): dev.digest(plan_path),
              str(Path(__file__)): dev.digest(Path(__file__))}
    for split, count in (("train", 16), ("calibration", 4), ("exposed", 8)):
        split_data[split], split_paths[split] = [], []
        for i in range(count):
            path = root / ("results_v3" if split == "exposed" else "results_v5") / f"{'confirmation' if split == 'exposed' else split}_{i:02d}"
            data = list(dev.dataset(path / "result.json"))
            if split == "exposed":
                phase_path = root / "phase_ablation_v5" / f"exposed_cohort_{i:02d}.npz"
                expected = dev.load(root / "phase_development_ablation_v5.json")["cohorts"][i]
                dev.require_hashes({str(phase_path): expected["phase_features_sha256"]})
                data[3] = np.concatenate([data[3], np.load(phase_path)["features"]], axis=1)
                bounds[str(phase_path)] = dev.digest(phase_path)
            split_data[split].append(data)
            split_paths[split].append(path)
            bounds[str(path / "result.json")] = dev.digest(path / "result.json")
    fold_models = {}
    for fold in range(4):
        fold_models[fold] = train_model([d for i, d in enumerate(split_data["train"]) if i % 4 != fold])
        dev.joblib.dump(fold_models[fold], out / f"sibling_fold_{fold}.joblib")
    map_path = root.parents[1] / "external/dog-recombination/campbell2016_average.plink.map"
    maps = dev.read_map(map_path)
    augmented, baseline_predictions, contexts = {}, {}, {}
    for split, data in split_data.items():
        augmented[split], baseline_predictions[split], contexts[split] = [], [], []
        for index, (d, path) in enumerate(zip(data, split_paths[split], strict=True)):
            base = dev.predict(d[2], d[3], baseline_model, candidate["temperature"], candidate["sibling_pooling"])
            group_model = fold_models[index % 4] if split == "train" else baseline_model
            grouping = dev.predict(d[2], d[3], group_model, candidate["temperature"], False)
            samples = {sample for r in d[2] for sample in dev.key(r)}
            records, _ = dev.read_haplotype_segments(path / "segments.ibd.gz", samples)
            masks = {t: dev.read_hbd_masks(path / "segments.hbd.gz", t) for t in PARAMETERS["hbd_cutoffs_cm"]}
            features, audit = context_features(d[2], probabilities(grouping), list(dev.LABELS).index("siblings"), records, masks, maps)
            # Truth is consulted only after extraction, for a separate grouping audit.
            membership = {sample: frozenset(group) for group in audit["predicted_groups"] for sample in group}
            actual = {dev.key(r): r["relation"] for r in d[1]["pairs"]}
            audit["grouping_truth_audit"] = {"false_grouped_pairs": sum(membership[a] == membership[b] and label != "siblings" for (a, b), label in actual.items()),
                "missed_sibling_pairs": sum(membership[a] != membership[b] and label == "siblings" for (a, b), label in actual.items()),
                "out_of_fold_grouping": split == "train"}
            np.savez_compressed(out / f"context_{split}_{index:02d}.npz", features=features)
            dev.dump(out / f"groups_{split}_{index:02d}.json", audit)
            augmented[split].append([*d[:3], np.concatenate([d[3], features], axis=1), *d[4:]])
            baseline_predictions[split].append(base)
            contexts[split].append(features)
    model = train_model(augmented["train"])
    dev.joblib.dump(model, out / "diagnostic_model.joblib")
    cx = np.concatenate([d[3][d[4]] for d in augmented["calibration"]])
    cy = np.concatenate([d[5] for d in augmented["calibration"]])
    losses = {str(t): float(dev.log_loss(cy, dev.temperature_scale(model.predict_proba(cx), t), labels=model.classes_)) for t in PARAMETERS["temperature_grid"]}
    temperature = float(min(losses, key=losses.get))
    reports = {}
    for split, data in augmented.items():
        rows = [dev.predict(d[2], d[3], model, temperature, False) for d in data]
        if split == "exposed":
            expected = dev.load(root / "phase_development_ablation_v5.json")["phase"]
            for key, value in aggregate(data, baseline_predictions[split]).items():
                if key in {"accuracy", "balanced_accuracy", "ece_10_equal_width"}:
                    assert abs(value - expected[key]) < 1e-12, key
                else:
                    assert value == expected[key], key
        cohort_rows = [{"cohort_id": f"{split}_{i:02d}", "truth": d[1]["pairs"], "predictions": p}
                       for i, (d, p) in enumerate(zip(data, rows, strict=True))]
        dev.dump(out / f"scoring_{split}.json", {"labels": list(dev.LABELS), "cohorts": cohort_rows})
        stratified = defaultdict(lambda: {"truth": {"pairs": []}, "base": [], "context": []})
        flips = {"base_correct_to_wrong": 0, "base_wrong_to_correct": 0, "strong_base_correct_to_wrong": 0}
        for i, (d, base, new) in enumerate(zip(data, baseline_predictions[split], rows, strict=True)):
            actual = {dev.key(r): r["relation"] for r in d[1]["pairs"]}
            for j, (b, n) in enumerate(zip(base, new, strict=True)):
                label = actual.get(dev.key(b))
                if label is None:
                    continue
                flips["base_correct_to_wrong"] += int(b["relation"] == label and n["relation"] != label)
                flips["base_wrong_to_correct"] += int(b["relation"] != label and n["relation"] == label)
                flips["strong_base_correct_to_wrong"] += int(b["support"] >= .9 and b["relation"] == label and n["relation"] != label)
                hbd = d[3][j, 130]
                stratum = "hbd_under_005" if hbd < .05 else "hbd_005_to_015" if hbd < .15 else "hbd_at_least_015"
                for group in (stratum, "complex" if d[0]["cohort"]["complex_topology"] else "simple"):
                    for original, target in (({**b, "relation": label}, stratified[group]["truth"]["pairs"]),
                                              (b, stratified[group]["base"]), (n, stratified[group]["context"])):
                        target.append({**original, "indID1": f"c{i}_{original['indID1']}", "indID2": f"c{i}_{original['indID2']}"})
        reports[split] = {"base": aggregate(data, baseline_predictions[split]), "context": aggregate(data, rows), "flips": flips,
            "strata": {group: {"base": dev.score(values["truth"], values["base"]), "context": dev.score(values["truth"], values["context"])} for group, values in stratified.items()},
            "context_medians_by_label": {label: np.median(np.concatenate([x[[j for j, r in enumerate(d[2]) if dev.key(r) in {dev.key(t) for t in d[1]["pairs"] if t["relation"] == label}]] for x, d in zip(contexts[split], data, strict=True)]), axis=0).tolist() for label in dev.LABELS}}
    dev.require_hashes(bounds)
    report = {"parameters": PARAMETERS, "temperature": temperature, "calibration_losses": losses, "reports": reports,
        "inputs_sha256": bounds, "outputs_sha256": {p.name: dev.digest(p) for p in sorted(out.iterdir())},
        "limitations": ["All cohorts use sibling-rich finite generator templates; general pedigree coverage not established.",
                        "Exposed8 are development checks, not independent confirmation.", "Training pair model scores are in-sample; only sibling grouping is out-of-fold.",
                        "No final candidate or qualification created; context model may damage strong correct pair calls."], "submission_ready": False}
    dev.dump(out / "report.json", report)
    print(json.dumps({split: {k: v for k, v in value.items() if k in ("base", "context", "flips")} for split, value in reports.items()}, indent=2))


if __name__ == "__main__":
    main(ROOT)
