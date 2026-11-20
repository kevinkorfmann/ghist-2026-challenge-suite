"""Prospective independent whole-family Task 11 development, never final inference.

Only synthetic founder cohorts are accepted. The held-out confirmation split is
inaccessible to the run command until the development model has been frozen.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.metadata
import json
import os
import platform
import shutil
import subprocess
import sys
import traceback
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.cluster import AgglomerativeClustering
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.metrics import adjusted_rand_score, log_loss
from sklearn.utils.class_weight import compute_sample_weight

from scripts.cluster_task11_families import affinity_matrix, canonical_labels, cluster
from scripts.fit_task11_dog_relationships import (
    LABELS,
    feature_names,
    merge_intervals,
    read_hbd_masks,
    read_observed,
    subtract_intervals,
)
from scripts.prepare_hapibd_vcf import prepare_vcf
from scripts.reconstruct_bonsai_pedigree import (
    coverage_intervals,
    read_haplotype_segments,
    read_map,
)
from scripts.validate_task11_complete_pipeline import gamete, partition_metrics

ROOT = Path(__file__).resolve().parents[1]
VERSION = "task11-detected-ibd-pair-first-development-v4"
HBD = (2.0, 5.0, 10.0)
PHASE_MIN_CM = 10.0
COUNTS = {"train": 16, "calibration": 4, "confirmation": 8}
SEED_BASE = {"train": 2026120100, "calibration": 2026120200, "confirmation": 2026120400}
GATES = {
    "accuracy_min": 0.80, "balanced_accuracy_min": 0.75,
    "each_label_recall_min": 0.60, "bootstrap_balanced_lower95_min": 0.70,
    "exact_family_recovery_min": 0.95, "expected_calibration_error_max": 0.10,
    "missing_samples_max": 0, "missing_true_pairs_max": 0,
    "extra_cross_family_pairs_max": 0,
}
DEPTHS = {
    "siblings": (1, 1, 2), "half_sibling": (1, 1, 1),
    "avuncular": (1, 2, 2), "half_avuncular": (1, 2, 1),
    "avuncular_once_removed": (1, 3, 2),
    "avuncular_twice_removed": (1, 4, 2),
    "full_1st_cousin": (2, 2, 2), "half_1st_cousin": (2, 2, 1),
}


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def dump(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")


def load(path: Path):
    return json.loads(path.read_text())


def require_hashes(entries: dict[str, str]) -> None:
    for path, expected in entries.items():
        if digest(Path(path)) != expected:
            raise ValueError(f"Input/source changed: {path}")


def source_identity() -> dict[str, str]:
    files = sorted((ROOT / "scripts").glob("*.py")) + sorted((ROOT / "src").rglob("*.py"))
    return {str(p.relative_to(ROOT)): digest(p) for p in files}


def runtime_versions() -> dict:
    return {"python": platform.python_version(), **{
        package: importlib.metadata.version(package)
        for package in ("numpy", "scipy", "scikit-learn", "msprime", "tskit", "pandas", "joblib")
    }}


def imported_source_tree() -> dict:
    found = {}
    for name, module in list(sys.modules.items()):
        path = getattr(module, "__file__", None)
        if path and name.startswith(("scripts.", "ghist_")):
            path = Path(path).resolve()
            if not path.is_relative_to(ROOT):
                raise ValueError(f"Imported source outside immutable snapshot: {name}: {path}")
            found[name] = {"path": str(path.relative_to(ROOT)), "sha256": digest(path)}
    return found


def require_runtime(protocol: dict) -> None:
    imported_source_tree()
    if runtime_versions() != protocol["runtime_versions"]:
        raise ValueError("Runtime differs from frozen protocol")


def validate_protocol(protocol: dict) -> dict:
    if protocol.get("schema") != "task11-prospective-development-v4" or protocol.get("method_version") != VERSION:
        raise ValueError("Unknown protocol schema or method")
    if protocol.get("gates") != GATES or set(protocol.get("gates", {})) != set(GATES):
        raise ValueError("Confirmation gates differ from preregistered implementation")
    return protocol["gates"]


def make_protocol(out: Path, map_path: Path, jar: Path, clustering_policy: str = "spectral", selection_evidence: Path | None = None) -> dict:
    if out.exists():
        raise FileExistsError(out)
    cohorts = []
    for split, count in COUNTS.items():
        for index in range(count):
            cohorts.append({
                "split": split, "index": index, "seed": SEED_BASE[split] + index,
                "population_size": [500, 1000, 2000][index % 3],
                "missing_rate": [0.0, 0.005, 0.02][(index // 3) % 3],
                "min_branch_size": 2, "max_branch_size": 4,
                "complex_topology": split == "confirmation" and index >= 4,
            })
    if clustering_policy != "spectral" and selection_evidence is None:
        raise ValueError("A changed clustering method requires bound training-selection evidence")
    inputs = {str(map_path.resolve()): digest(map_path), str(jar.resolve()): digest(jar)}
    if selection_evidence is not None:
        inputs[str(selection_evidence.resolve())] = digest(selection_evidence)
    payload = {
        "schema": "task11-prospective-development-v4",
        "method_version": VERSION, "cohorts": cohorts, "gates": GATES,
        "clustering_policy": clustering_policy,
        "clustering_selection": "Select spectral, single-link or average-link using 16 training cohorts only, by exact family recovery followed by fewer missing+extra true pairs. Fixed 10cM segment threshold, supplied 8 families; 20cM is diagnostic sensitivity only. Training tie between single/average resolves to average to reduce chaining.",
        "inputs": inputs,
        "clustering_selection_evidence": str(selection_evidence.resolve()) if selection_evidence else None,
        "source_sha256": source_identity(), "labels": list(LABELS), "runtime_versions": runtime_versions(),
        "imported_sources": imported_source_tree(),
        "detector": {"min_seed": 1.0, "min_extend": 0.2, "min_output": 2.0,
                     "min_markers": 50, "max_gap": 1000},
        "founders": "Independent whole-cohort msprime coalescent samples from one population; Ne 500/1000/2000; species recombination map; mutation 1.25e-8; biallelic sites with minor allele count >=2 and spacing >=500bp. No empirical breed distribution assumed.",
        "pedigrees": "8 unrelated focal families, one template per accepted label; two 2-4 member sibships plus optional extra observed siblings. Confirmation cohorts4-7 include new third branches spanning avuncular/once-removed/twice-removed labels; all relatives gene-dropped jointly. Inbreeding above focal shared parent couple cycles unrelated/first-cousin/full-sibling mating; direct focal relationship is truth even with background ancestral loops.",
        "preprocessing": "Known true phase; MCAR complete-genotype masking; drop any site with missing genotype; no imputation.",
        "model": "ExtraTrees 800 trees, leaf3, max_features0.75, class-balanced weights; concatenated 2/5/10cM HBD feature vectors plus endpoint HBD fractions; whole-cohort fitting split.",
        "phase_extension": "One fixed >=10cM phase variant after the >=2cM extension failed development recall. Union each homolog edge's duplicate/overlapping intervals before the >=10cM length filter. Genuine IBD2 (two-edge matching) is masked using ALL detected segments, including shorter ones. At HBD cutoffs2/5/10cM mask union of both samples' HBD. For each endpoint sum chromosome-wise larger homolog mass / total retained homolog mass. Sort endpoint scores, retained cM denominators, and retained/pre-HBD denominators into18 features. Zero exposure has score0, denominator0, fraction0. Chromosome homolog orientation and pair endpoint order are arbitrary; features must be invariant. No sex, age, labels, family truth or sample names enter features. All original detector/clustering/classifier hyperparameters unchanged. Matched no-phase ablation uses identical development cohorts. Exposed seeds2026120300-307 are development-only and never training or new confirmation. Motivation for10cM is independent same-genome exact/detected Task10 evidence of short opposite-homolog contamination; transfer to dogs is a hypothesis, not a fitted threshold search.",
        "calibration": "Choose temperature from [0.75,1,1.25,1.5,2] by calibration-cohort multiclass log loss; no confirmation selection.",
        "global_policy": "Eligible complete-link sibling components require every internal sibling probability>=0.95. Pool probabilities across equivalent full-sibling component edges only for original confidence<0.90 and pooled confidence>=0.90. Enable only if calibration balanced accuracy improves and introduces no errors in previously correct calls. Never Bonsai override.",
        "confirmation_policy": "Freeze model, temperature and global policy before any confirmation execution; 8 entirely new cohorts; unconditional truth-pair scoring counts omissions wrong and reports cross-family extras. Cohort-bootstrap 2000 draws, seed20301120.",
        "confidence": "Temperature-scaled model support, evaluated by fixed 10 equal-width ECE bins; not assumed calibrated before confirmation.",
        "submission_ready": False, "final_challenge_access": False, "submission_count": 0,
    }
    dump(out, payload)
    return payload


def focal_pedigree(label: str, sizes: tuple[int, int], third_branch: bool = False, complex_topology: bool = False) -> dict:
    """Build observed terminal relatives with explicit, non-genotype-derived truth."""
    left_depth, right_depth, shared = DEPTHS[label]
    parents: dict[str, tuple[str, str] | None] = {}
    counter = 0

    def node(pair=None):
        nonlocal counter
        name = f"n{counter}"
        counter += 1
        parents[name] = pair
        return name

    root1, root2 = node(), node()
    other = root2 if shared == 2 else node()

    def branch(depth, pair, size):
        for _ in range(depth - 1):
            pair = (node(pair), node())
        return [node(pair) for _ in range(size)]

    left = branch(left_depth, (root1, root2), sizes[0])
    right = branch(right_depth, (root1, other), sizes[1])
    groups = [left, right]
    rows = []
    for group in groups:
        rows.extend((a, b, "siblings") for a, b in combinations(group, 2))
    rows.extend((a, b, label) for a in left for b in right)
    if third_branch:
        # New confirmation topologies put relatives at several depths in one
        # observed family while keeping every true pair in the accepted labels.
        cross_left, cross_right = label, "siblings"
        if complex_topology and label in {"avuncular", "avuncular_once_removed"}:
            hidden = node(parents[right[0]])
            mate = node()
            extra = [node((hidden, mate)) for _ in range(2)]
            cross_left = {"avuncular": "avuncular_once_removed", "avuncular_once_removed": "avuncular_twice_removed"}[label]
            cross_right = "avuncular"
        elif complex_topology and label == "avuncular_twice_removed":
            hidden = parents[right[0]][0]
            extra = [node(parents[hidden]) for _ in range(2)]
            cross_left, cross_right = "avuncular_once_removed", "avuncular"
        else:
            extra = [node(parents[right[0]]) for _ in range(2)]
        rows.extend((a, b, cross_right) for a in right for b in extra)
        rows.extend((a, b, "siblings") for a, b in combinations(extra, 2))
        rows.extend((a, b, cross_left) for a in left for b in extra)
        groups.append(extra)
    observed = [n for group in groups for n in group]
    return {"parents": parents, "observed": observed, "pairs": rows,
            "central_founders": [root1, root2], "label": label}


def expand_inbreeding(pedigree: dict, kind: str) -> dict:
    """Add background ancestry without changing focal relationship truth."""
    original = pedigree["parents"]
    prefix = {}
    if kind == "unrelated":
        return dict(original)
    a, b = pedigree["central_founders"]
    prefix["bg0"] = prefix["bg1"] = None
    if kind == "siblings":
        prefix[a] = prefix[b] = ("bg0", "bg1")
    elif kind == "first_cousins":
        prefix["bg2"] = prefix["bg3"] = None
        prefix["bg4"] = prefix["bg5"] = ("bg0", "bg1")
        prefix[a] = ("bg4", "bg2")
        prefix[b] = ("bg5", "bg3")
    else:
        raise ValueError(kind)
    for name, pair in original.items():
        if name not in prefix:
            prefix[name] = pair
    return prefix


def generate_genotypes(map_path: Path, out: Path, config: dict) -> dict:
    import msprime

    rng = np.random.default_rng(config["seed"])
    maps = read_map(map_path)
    records = []
    samples, families, truth_pairs = [], [], []
    next_sample = 0
    kinds = ("unrelated", "first_cousins", "siblings")
    for index, label in enumerate(LABELS):
        sizes = tuple(int(x) for x in rng.integers(config["min_branch_size"], config["max_branch_size"] + 1, 2))
        complex_topology = bool(config.get("complex_topology", False))
        third = bool((config["index"] + index) % 2) or (complex_topology and label.startswith("avuncular"))
        pedigree = focal_pedigree(label, sizes, third, complex_topology)
        parents = expand_inbreeding(pedigree, kinds[(index + config["index"]) % 3])
        names = {n: f"synthetic11_{next_sample + i}" for i, n in enumerate(pedigree["observed"])}
        next_sample += len(names)
        families.append(list(names.values()))
        samples.extend(names.values())
        truth_pairs.extend({"indID1": names[a], "indID2": names[b], "relation": relation,
                            "family": index} for a, b, relation in pedigree["pairs"])
        records.append({"parents": parents, "observed": pedigree["observed"],
                        "kind": kinds[(index + config["index"]) % 3], "names": names})
    founder_count = sum(sum(pair is None for pair in r["parents"].values()) for r in records)
    hbd_cm = np.zeros(len(samples))
    genome_cm, markers, missing_calls = 0.0, 0, 0
    with gzip.open(out / "synthetic.vcf.gz", "wt") as destination:
        destination.write("##fileformat=VCFv4.2\n")
        destination.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Known simulated phase">\n')
        for chrom, (bp, _) in maps.items():
            destination.write(f"##contig=<ID={chrom},length={int(bp[-1]) + 1}>\n")
        destination.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
        for chrom, (map_bp, map_cm) in maps.items():
            positions = np.concatenate(([0.0], map_bp.astype(float), [float(map_bp[-1] + 1)]))
            # Remove duplicate 0 if the original map begins there.
            positions = np.unique(positions)
            cm_boundaries = np.interp(positions, map_bp, map_cm)
            rates = np.diff(cm_boundaries) / 100.0 / np.diff(positions)
            ts = msprime.sim_ancestry(samples=founder_count, population_size=config["population_size"],
                recombination_rate=msprime.RateMap(position=positions, rate=rates),
                random_seed=int(rng.integers(1, 2**31)))
            ts = msprime.sim_mutations(ts, rate=1.25e-8, random_seed=int(rng.integers(1, 2**31)))
            site_bp, alleles = [], []
            last = -500
            for variant in ts.variants():
                bp = int(variant.site.position)
                g = variant.genotypes
                if (bp < map_bp[0] or bp > map_bp[-1] or bp - last < 500
                    or len(variant.alleles) != 2 or min(int(g.sum()), len(g) - int(g.sum())) < 2):
                    continue
                site_bp.append(bp)
                alleles.append(g.copy())
                last = bp
            if len(site_bp) < 100:
                raise ValueError(f"Insufficient segregating markers on chromosome {chrom}")
            bp = np.array(site_bp)
            cm = np.interp(bp, map_bp, map_cm)
            founder_alleles = np.asarray(alleles, dtype=np.int8).T
            offset, ancestry_samples = 0, []
            for record in records:
                ancestry = {}
                for name, pair in record["parents"].items():
                    if pair is None:
                        ancestry[name] = np.repeat(np.arange(offset, offset + 2)[:, None], len(bp), axis=1)
                        offset += 2
                    else:
                        ancestry[name] = np.stack([gamete(ancestry[p], cm, rng) for p in pair])
                ancestry_samples.extend(ancestry[n] for n in record["observed"])
            ancestors = np.stack(ancestry_samples)
            genotypes = founder_alleles[ancestors, np.arange(len(bp))]
            is_hbd = ancestors[:, 0] == ancestors[:, 1]
            hbd_cm += (is_hbd[:, :-1] * np.diff(cm)).sum(axis=1)
            genome_cm += float(map_cm[-1] - map_cm[0])
            missing = rng.random((len(samples), len(bp))) < config["missing_rate"]
            missing_calls += int(missing.sum())
            for j, position in enumerate(bp):
                calls = [".|." if missing[i, j] else f"{int(g[0])}|{int(g[1])}"
                         for i, g in enumerate(genotypes[:, :, j])]
                destination.write(f"{chrom}\t{position}\t.\tA\tC\t.\tPASS\t.\tGT\t" + "\t".join(calls) + "\n")
            markers += len(bp)
    truth = {"samples": samples, "families": families, "pairs": truth_pairs,
             "pedigrees": records, "genome_cm": genome_cm, "markers": markers,
             "founder_individuals": founder_count, "seed": config["seed"],
             "missing_genotypes": missing_calls, "true_hbd_cm": dict(zip(samples, hbd_cm.tolist()))}
    dump(out / "truth.json", truth)
    return truth


def execute(name: str, command: list[str], out: Path, stages: list[dict]) -> None:
    with (out / f"{name}.out").open("w") as stdout, (out / f"{name}.err").open("w") as stderr:
        result = subprocess.run(command, cwd=ROOT, stdout=stdout, stderr=stderr, check=False)
    stages.append({"stage": name, "exit_code": result.returncode, "command": command})
    if result.returncode:
        raise RuntimeError(f"{name} failed with {result.returncode}")


def pair_haplotype_scores(records_by_chromosome: dict, maps: dict, masks: dict) -> np.ndarray:
    """IBD2/HBD-masked chromosome-wise homolog concentration, with exposure.

    The concentration follows PONDEROSA equation2; unioned coverage makes it
    invariant to redundant calls. Endpoint scores and exposures are sorted.
    Zero exposure is encoded as zero, with explicit zero denominator features.
    """
    numerator, denominator, before_hbd = np.zeros(2), np.zeros(2), np.zeros(2)
    for chromosome, records in records_by_chromosome.items():
        if chromosome not in maps:
            raise ValueError("Phase feature chromosome absent from map")
        bp, cm = maps[chromosome]
        if any(start < bp[0] or end > bp[-1] for start, end, _, _ in records):
            raise ValueError("Phase feature interval exceeds map bounds")
        full = [(a, b) for a, b, is_full in coverage_intervals(records) if is_full]
        excluded = merge_intervals(full + masks.get(chromosome, []))
        by_edge = defaultdict(list)
        for start, end, hap1, hap2 in records:
            by_edge[(hap1, hap2)].append((start, end))
        eligible = [(start, end, hap1, hap2) for (hap1, hap2), intervals in by_edge.items()
                    for start, end in merge_intervals(intervals)
                    if float(np.diff(np.interp([start, end], bp, cm))[0]) >= PHASE_MIN_CM]
        after, before = defaultdict(list), defaultdict(list)
        for start, end, hap1, hap2 in eligible:
            for side, haplotype in ((0, hap1), (1, hap2)):
                after[(side, haplotype)].extend(subtract_intervals(start, end, excluded))
                before[(side, haplotype)].extend(subtract_intervals(start, end, full))
        for side in (0, 1):
            masses = []
            for haplotype in (1, 2):
                def length(intervals, positions=bp, centimorgans=cm):
                    return sum(float(np.diff(np.interp([a, b], positions, centimorgans))[0])
                               for a, b in merge_intervals(intervals))
                masses.append(length(after[(side, haplotype)]))
                before_hbd[side] += length(before[(side, haplotype)])
            numerator[side] += max(masses)
            denominator[side] += sum(masses)
    scores = np.divide(numerator, denominator, out=np.zeros(2), where=denominator > 0)
    fractions = np.divide(denominator, before_hbd, out=np.zeros(2), where=before_hbd > 0)
    return np.array([*sorted(scores), *sorted(denominator), *sorted(fractions)])


def haplotype_features(ibd_path: Path, hbd_path: Path, map_path: Path, pairs: list) -> np.ndarray:
    grouped, _ = read_haplotype_segments(ibd_path, {sample for pair in pairs for sample in pair})
    maps = read_map(map_path)
    by_pair = defaultdict(dict)
    for (a, b, chromosome), records in grouped.items():
        by_pair[(a, b)][chromosome] = records
    matrices = []
    for threshold in HBD:
        all_masks = read_hbd_masks(hbd_path, threshold)
        rows = []
        for a, b in pairs:
            masks = {chromosome: merge_intervals(all_masks.get((a, chromosome), []) +
                     all_masks.get((b, chromosome), [])) for chromosome in maps}
            rows.append(pair_haplotype_scores(by_pair[tuple(sorted((a, b)))], maps, masks))
        matrices.append(np.asarray(rows))
    return np.concatenate(matrices, axis=1)


def phase_feature_names() -> list[str]:
    return [f"hbd{int(t)}_phase_{name}" for t in HBD for name in
            ("concentration_low", "concentration_high", "retained_cm_low", "retained_cm_high",
             "retained_fraction_low", "retained_fraction_high")]


def detect_cluster_extract(out: Path, map_path: Path, jar: Path, samples: list[str], policy: str) -> dict:
    stages = []
    preparation = prepare_vcf(out / "synthetic.vcf.gz", out / "prepared.vcf.gz", missing_policy="drop-site")
    dump(out / "preparation.json", preparation)
    execute("hapibd", ["java", "-Xmx12g", "-jar", str(jar), f"gt={out / 'prepared.vcf.gz'}",
        f"map={map_path}", f"out={out / 'segments'}", "min-seed=1.0", "min-extend=0.2",
        "min-output=2.0", "min-markers=50", "max-gap=1000", "nthreads=4"], out, stages)
    execute("summary", [sys.executable, str(ROOT / "scripts/summarize_hapibd.py"),
        "--ibd", str(out / "segments.ibd.gz"), "--hbd", str(out / "segments.hbd.gz"),
        "--map", str(map_path), "--out-dir", str(out / "summary")], out, stages)
    frame = pd.read_csv(out / "summary/pair_features.csv")
    solutions = []
    for threshold in (10, 20):
        affinity = affinity_matrix(frame, samples, f"relatedness_ge_{threshold}")
        if policy == "spectral":
            solutions.extend(canonical_labels(cluster(affinity, 8, seed), samples) for seed in range(8))
        elif policy in {"single", "average"}:
            distance = affinity.max() - affinity
            np.fill_diagonal(distance, 0)
            labels = AgglomerativeClustering(n_clusters=8, metric="precomputed", linkage=policy).fit_predict(distance)
            solutions.append(canonical_labels(labels, samples))
        else:
            raise ValueError("Unknown clustering policy")
    reference = solutions[0]
    families = [[s for s, assigned in zip(samples, reference, strict=True) if assigned == i] for i in range(8)]
    clusters = {"families": families, "minimum_adjusted_rand": min(float(adjusted_rand_score(reference, s)) for s in solutions),
                "family_count_supplied": 8, "reference": "ge10", "policy": policy, "no_truth_partition_used": True}
    dump(out / "clusters.json", clusters)
    matrices, reference_pairs, reference_indices = [], None, None
    for threshold in HBD:
        matrix, pairs, indices = read_observed(out / "segments.ibd.gz", out / "segments.hbd.gz", map_path,
                                             out / "clusters.json", threshold)
        if reference_pairs is not None and pairs != reference_pairs:
            raise ValueError("HBD feature rows differ")
        reference_pairs, reference_indices = pairs, indices
        matrices.append(matrix)
    hbd = pd.read_csv(out / "summary/hbd_features.csv").set_index("individual")
    fractions = {sample: float(hbd.loc[sample, "inbreeding_fraction"]) if sample in hbd.index else 0.0 for sample in samples}
    endpoint = np.array([sorted((fractions[a], fractions[b])) for a, b in reference_pairs])
    phase = haplotype_features(out / "segments.ibd.gz", out / "segments.hbd.gz", map_path, reference_pairs)
    features = np.concatenate([*matrices, endpoint, phase], axis=1)
    if not np.isfinite(features).all():
        raise ValueError("Nonfinite features")
    np.savez_compressed(out / "features.npz", features=features)
    dump(out / "feature_pairs.json", [{"indID1": a, "indID2": b, "family": int(family)}
         for (a, b), family in zip(reference_pairs, reference_indices, strict=True)])
    return {"stages": stages, "clusters": clusters, "features": len(features),
            "phase_exposure": {str(t): {"zero_endpoint_denominators": int((phase[:, i * 6 + 2:i * 6 + 4] == 0).sum()),
                               "retained_cm_quantiles": np.quantile(phase[:, i * 6 + 2:i * 6 + 4], [0, .25, .5, .75, 1]).tolist(),
                               "retained_fraction_quantiles": np.quantile(phase[:, i * 6 + 4:i * 6 + 6], [0, .25, .5, .75, 1]).tolist()}
                               for i, t in enumerate(HBD)},
            "feature_names": [f"hbd{int(t)}_{n}" for t in HBD for n in feature_names()] + ["min_endpoint_hbd", "max_endpoint_hbd"] + phase_feature_names()}


def key(row: dict) -> tuple[str, str]:
    return tuple(sorted((row["indID1"], row["indID2"])))


def score(truth: dict, rows: list[dict]) -> dict:
    expected = {key(r): r for r in truth["pairs"]}
    predicted = {}
    for row in rows:
        if key(row) in predicted or row["indID1"] == row["indID2"]:
            raise ValueError("Duplicate or self-pair")
        predicted[key(row)] = row
    counts, correct, confusion = Counter(), Counter(), defaultdict(Counter)
    confidence, outcomes = [], []
    for pair, actual in expected.items():
        label = actual["relation"]
        counts[label] += 1
        p = predicted.get(pair, {})
        observed = p.get("relation", "MISSING")
        correct[label] += observed == label
        confusion[label][observed] += 1
        if pair in predicted:
            confidence.append(float(p["support"]))
            outcomes.append(observed == label)
    recall = {label: correct[label] / counts[label] if counts[label] else 0.0 for label in LABELS}
    ece = 0.0
    bins = np.minimum((np.asarray(confidence) * 10).astype(int), 9)
    for bin_index in range(10):
        selected = np.flatnonzero(bins == bin_index).tolist()
        if selected:
            ece += len(selected) / len(confidence) * abs(np.mean([confidence[i] for i in selected]) - np.mean([outcomes[i] for i in selected]))
    return {"truth_pairs": len(expected), "correct": sum(correct.values()),
            "accuracy": sum(correct.values()) / len(expected), "balanced_accuracy": float(np.mean(list(recall.values()))),
            "recall": recall, "counts": dict(counts), "correct_counts": dict(correct),
            "confusion": {label: dict(count) for label, count in confusion.items()},
            "missing_true_pairs": len(expected.keys() - predicted.keys()),
            "extra_cross_family_pairs": len(predicted.keys() - expected.keys()),
            "ece_10_equal_width": float(ece)}


def temperature_scale(probabilities: np.ndarray, temperature: float) -> np.ndarray:
    values = np.maximum(probabilities, 1e-12) ** (1.0 / temperature)
    return values / values.sum(axis=1, keepdims=True)


def sibling_pool(rows: list[dict], probabilities: np.ndarray) -> np.ndarray:
    """Apply only complete-link full-sibling equivalence, never override strong calls."""
    result = probabilities.copy()
    lookup = {key(r): i for i, r in enumerate(rows)}
    families = defaultdict(set)
    for row in rows:
        families[row["family"]].update(key(row))
    sibling_index = list(LABELS).index("siblings")
    for members in families.values():
        components = [{s} for s in sorted(members)]
        changed = True
        while changed:
            changed = False
            for i, a in enumerate(components):
                for j in range(i + 1, len(components)):
                    b = components[j]
                    if all(probabilities[lookup[tuple(sorted((x, y)))], sibling_index] >= 0.95 for x in a for y in b):
                        components[i] = a | b
                        components.pop(j)
                        changed = True
                        break
                if changed:
                    break
        for a, b in combinations(components, 2):
            indices = [lookup[tuple(sorted((x, y)))] for x in a for y in b]
            if len(indices) < 2:
                continue
            pooled = probabilities[indices].mean(axis=0)
            if pooled.max() >= 0.90:
                for index in indices:
                    if probabilities[index].max() < 0.90:
                        result[index] = pooled
    return result


def predict(rows: list[dict], features: np.ndarray, model, temperature: float, pooling: bool) -> list[dict]:
    raw = model.predict_proba(features)
    probabilities = np.zeros((len(raw), len(LABELS)))
    for index, label in enumerate(model.classes_):
        probabilities[:, list(LABELS).index(str(label))] = raw[:, index]
    probabilities = temperature_scale(probabilities, temperature)
    if pooling:
        probabilities = sibling_pool(rows, probabilities)
    return [{**row, "relation": LABELS[int(np.argmax(p))], "support": float(p.max()),
             "probabilities": dict(zip(LABELS, p.tolist()))} for row, p in zip(rows, probabilities, strict=True)]


def run_cohort(args) -> dict:
    protocol = load(args.protocol)
    validate_protocol(protocol)
    config = next(c for c in protocol["cohorts"] if c["split"] == args.split and c["index"] == args.index)
    if str(args.map.resolve()) not in protocol["inputs"] or str(args.hapibd.resolve()) not in protocol["inputs"]:
        raise ValueError("Unregistered map or detector")
    source_bound = {str(ROOT / path): value for path, value in protocol["source_sha256"].items()}
    bounds = {**protocol["inputs"], **source_bound, str(args.protocol.resolve()): digest(args.protocol)}
    if args.split == "confirmation":
        if args.baseline is not None:
            raise ValueError("Confirmation cannot replay a development cohort")
        if args.candidate is None:
            raise ValueError("Confirmation requires prospectively frozen candidate")
        candidate = load(args.candidate)
        if candidate["method_version"] != VERSION or not candidate["frozen_before_confirmation"]:
            raise ValueError("Unfrozen or unknown candidate")
        if str(Path(candidate["model_path"]).resolve()) not in candidate["bound_inputs"]:
            raise ValueError("Unbound candidate model")
        if candidate["protocol_sha256"] != digest(args.protocol):
            raise ValueError("Candidate protocol mismatch")
        bounds.update(candidate["bound_inputs"])
        bounds[str(args.candidate.resolve())] = digest(args.candidate)
    require_runtime(protocol)
    require_hashes(bounds)
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    dump(out / "protocol.json", {"master_sha256": digest(args.protocol), "cohort": config,
        "source_sha256": protocol["source_sha256"], "inputs": protocol["inputs"],
        "candidate_sha256": digest(args.candidate) if args.candidate else None,
        "job_id": os.environ.get("SLURM_JOB_ID"), "python": platform.python_version()})
    result = {"cohort": config, "protocol_sha256": digest(out / "protocol.json"),
              "submission_ready": False, "submission_count": 0}
    try:
        if args.baseline is None:
            truth = generate_genotypes(args.map, out, config)
        else:
            baseline_result = load(args.baseline / "result.json")
            baseline_config = {"complex_topology": False, **baseline_result["cohort"]}
            if not baseline_result.get("completed") or baseline_config != config:
                raise ValueError("Baseline failed or cohort identity changed")
            baseline_bound = {str(args.baseline / name): h for name, h in baseline_result["output_sha256"].items()}
            require_hashes(baseline_bound)
            bounds.update(baseline_bound)
            bounds[str(args.baseline / "result.json")] = digest(args.baseline / "result.json")
            for name in ("synthetic.vcf.gz", "truth.json"):
                shutil.copyfile(args.baseline / name, out / name)
            truth = load(out / "truth.json")
            result["development_genotype_replay"] = {"baseline": str(args.baseline), "result_sha256": digest(args.baseline / "result.json"), "reason": "Training-only family clustering comparison; original immutable genotype/truth retained, detector and downstream stages rerun."}
        bounds.update({str(out / name): digest(out / name) for name in ("synthetic.vcf.gz", "truth.json")})
        inference = detect_cluster_extract(out, args.map, args.hapibd, truth["samples"], protocol["clustering_policy"])
        result["inference"] = inference
        result["partition"] = partition_metrics(truth["families"], inference["clusters"]["families"])
        if args.split == "confirmation":
            model = joblib.load(Path(candidate["model_path"]))
            rows = predict(load(out / "feature_pairs.json"), np.load(out / "features.npz")["features"],
                model, candidate["temperature"], candidate["sibling_pooling"])
            dump(out / "predictions.json", rows)
            result["score"] = score(truth, rows)
        require_hashes(bounds)
        result["completed"] = True
    except Exception as error:  # noqa: BLE001 - retain complete failed-cohort audit; main exits nonzero
        result.update(completed=False, failure={"type": type(error).__name__, "message": str(error), "traceback": traceback.format_exc()})
    result["output_sha256"] = {str(path.relative_to(out)): digest(path) for path in sorted(out.rglob("*")) if path.is_file()}
    dump(out / "result.json", result)
    return result


def dataset(result_path: Path):
    result = load(result_path)
    if not result.get("completed"):
        raise ValueError(f"Incomplete cohort {result_path}")
    out = result_path.parent
    require_hashes({str(out / p): h for p, h in result["output_sha256"].items()})
    truth, rows = load(out / "truth.json"), load(out / "feature_pairs.json")
    lookup = {key(r): r["relation"] for r in truth["pairs"]}
    selected = [i for i, row in enumerate(rows) if key(row) in lookup]
    # False cross-family rows never become invented training relationships.
    features = np.load(out / "features.npz")["features"]
    return result, truth, rows, features, selected, np.array([lookup[key(rows[i])] for i in selected])


def fit(args) -> dict:
    protocol = load(args.protocol)
    validate_protocol(protocol)
    require_runtime(protocol)
    require_hashes({**protocol["inputs"], **{str(ROOT / p): h for p, h in protocol["source_sha256"].items()}})
    args.out.mkdir(parents=True, exist_ok=False)
    train = [dataset(args.results / f"train_{i:02d}/result.json") for i in range(COUNTS["train"])]
    calibration = [dataset(args.results / f"calibration_{i:02d}/result.json") for i in range(COUNTS["calibration"])]
    for split, cohorts in (("train", train), ("calibration", calibration)):
        for i in range(COUNTS[split]):
            if load(args.results / f"{split}_{i:02d}/protocol.json")["master_sha256"] != digest(args.protocol):
                raise ValueError("Development protocol mismatch")
        if [c[0]["cohort"]["seed"] for c in cohorts] != [SEED_BASE[split] + i for i in range(COUNTS[split])]:
            raise ValueError("Unexpected development seed identities")
    x = np.concatenate([d[3][d[4]] for d in train])
    y = np.concatenate([d[5] for d in train])
    model = ExtraTreesClassifier(n_estimators=800, min_samples_leaf=3, max_features=0.75,
                                n_jobs=4, random_state=20300111)
    model.fit(x, y, sample_weight=compute_sample_weight("balanced", y))
    if set(model.classes_) != set(LABELS):
        raise ValueError("Training omitted accepted labels")
    losses = {}
    cx = np.concatenate([d[3][d[4]] for d in calibration])
    cy = np.concatenate([d[5] for d in calibration])
    for temperature in (0.75, 1.0, 1.25, 1.5, 2.0):
        losses[str(temperature)] = float(log_loss(cy, temperature_scale(model.predict_proba(cx), temperature), labels=model.classes_))
    temperature = min(losses, key=losses.get)
    original, pooled, truths = [], [], {"pairs": []}
    old_correct_flips = 0
    for index, d in enumerate(calibration):
        raw = predict(d[2], d[3], model, float(temperature), False)
        global_rows = predict(d[2], d[3], model, float(temperature), True)
        actual = {key(r): r["relation"] for r in d[1]["pairs"]}
        old_correct_flips += sum(actual.get(key(a)) == a["relation"] and actual.get(key(b)) != b["relation"] for a, b in zip(raw, global_rows, strict=True))
        for source, dest in ((raw, original), (global_rows, pooled), (d[1]["pairs"], truths["pairs"])):
            dest.extend({**r, "indID1": f"c{index}_{r['indID1']}", "indID2": f"c{index}_{r['indID2']}"} for r in source)
    raw_score, pooled_score = score(truths, original), score(truths, pooled)
    enable = pooled_score["balanced_accuracy"] > raw_score["balanced_accuracy"] and old_correct_flips == 0
    model_path = args.out / "model.joblib"
    joblib.dump(model, model_path)
    bound = {str(args.protocol.resolve()): digest(args.protocol), str(model_path.resolve()): digest(model_path)}
    for split in ("train", "calibration"):
        for i in range(COUNTS[split]):
            p = args.results / f"{split}_{i:02d}/result.json"
            bound[str(p.resolve())] = digest(p)
    candidate = {"method_version": VERSION, "protocol_sha256": digest(args.protocol),
                 "model_path": str(model_path.resolve()), "temperature": float(temperature),
                 "sibling_pooling": enable, "calibration_log_losses": losses,
                 "calibration_pair_score": raw_score, "calibration_pool_score": pooled_score,
                 "pooling_new_errors_on_previously_correct": int(old_correct_flips),
                 "bound_inputs": bound, "training_pairs": len(y), "training_label_counts": dict(Counter(y)),
                 "frozen_before_confirmation": True, "submission_ready": False,
                 "source_sha256": protocol["source_sha256"], "gates": protocol["gates"]}
    dump(args.out / "candidate.json", candidate)
    return candidate


def summarize(args) -> dict:
    protocol = load(args.protocol)
    gates = validate_protocol(protocol)
    datasets = [dataset(args.results / f"confirmation_{i:02d}/result.json") for i in range(COUNTS["confirmation"])]
    candidate_hashes = set()
    for i, d in enumerate(datasets):
        p = load(args.results / f"confirmation_{i:02d}/protocol.json")
        if p["master_sha256"] != digest(args.protocol) or d[0]["cohort"]["seed"] != SEED_BASE["confirmation"] + i:
            raise ValueError("Confirmation protocol or seed mismatch")
        candidate_hashes.add(p["candidate_sha256"])
    if len(candidate_hashes) != 1 or None in candidate_hashes:
        raise ValueError("Confirmation must use one frozen candidate")
    crosscheck = {"labels": list(LABELS), "cohorts": [
        {"cohort_id": str(d[0]["cohort"]["seed"]), "truth": d[1]["pairs"],
         "predictions": load(args.results / f"confirmation_{i:02d}/predictions.json")}
        for i, d in enumerate(datasets)
    ]}
    dump(args.out.with_name("independent_scoring_input.json"), crosscheck)
    truth, rows = {"pairs": []}, []
    for index, d in enumerate(datasets):
        for source, dest in ((d[1]["pairs"], truth["pairs"]), (load(args.results / f"confirmation_{index:02d}/predictions.json"), rows)):
            dest.extend({**r, "indID1": f"c{index}_{r['indID1']}", "indID2": f"c{index}_{r['indID2']}"} for r in source)
    overall = score(truth, rows)
    recomputed_scores = [score(d[1], load(args.results / f"confirmation_{i:02d}/predictions.json")) for i, d in enumerate(datasets)]
    recomputed_partitions = [partition_metrics(d[1]["families"], load(args.results / f"confirmation_{i:02d}/clusters.json")["families"]) for i, d in enumerate(datasets)]
    for d, measured in zip(datasets, recomputed_scores, strict=True):
        if measured != d[0]["score"]:
            raise ValueError("Stored confirmation score disagrees with retained truth/predictions")
    rng = np.random.default_rng(20301120)
    recalls = []
    for _ in range(2000):
        selected = rng.integers(0, len(datasets), len(datasets))
        correct, counts = Counter(), Counter()
        for i in selected:
            correct.update(recomputed_scores[i]["correct_counts"])
            counts.update(recomputed_scores[i]["counts"])
        recalls.append(np.mean([correct[label] / counts[label] for label in LABELS]))
    overall["cohort_bootstrap_balanced_interval95"] = np.quantile(recalls, [0.025, 0.975]).tolist()
    exact = sum(p["exactly_recovered_families"] for p in recomputed_partitions) / (8 * len(datasets))
    missing = sum(len(p["missing_samples"]) for p in recomputed_partitions)
    checks = {"accuracy": overall["accuracy"] >= gates["accuracy_min"],
        "balanced_accuracy": overall["balanced_accuracy"] >= gates["balanced_accuracy_min"],
        "each_label_recall": min(overall["recall"].values()) >= gates["each_label_recall_min"],
        "bootstrap_lower": overall["cohort_bootstrap_balanced_interval95"][0] >= gates["bootstrap_balanced_lower95_min"],
        "family_recovery": exact >= gates["exact_family_recovery_min"],
        "missing_samples": missing == 0, "missing_true_pairs": overall["missing_true_pairs"] == 0,
        "extra_pairs": overall["extra_cross_family_pairs"] == 0,
        "calibration": overall["ece_10_equal_width"] <= gates["expected_calibration_error_max"]}
    summary = {"method_version": VERSION, "score": overall, "checks": checks,
               "protocol_sha256": digest(args.protocol), "candidate_sha256": next(iter(candidate_hashes)),
               "confirmation_results_sha256": {str(args.results / f"confirmation_{i:02d}/result.json"): digest(args.results / f"confirmation_{i:02d}/result.json") for i in range(len(datasets))},
               "crosscheck_input_sha256": digest(args.out.with_name("independent_scoring_input.json")),
               "predeclared_criteria_passed": all(checks.values()), "exact_family_recovery": exact,
               "cohorts": [d[0]["cohort"] for d in datasets], "submission_ready": False,
               "submission_count": 0, "limitations": ["No final challenge inference performed.", "Four confirmation cohorts contain new three-branch generation structures; finite templates still do not cover every possible family topology.", "Coalescent population/missingness stresses cover prespecified assumptions, not all possible dog histories."]}
    dump(args.out, summary)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("protocol")
    p.add_argument("--map", type=Path, required=True)
    p.add_argument("--hapibd", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--clustering-policy", choices=("spectral", "single", "average"), required=True)
    p.add_argument("--selection-evidence", type=Path)
    p = sub.add_parser("run")
    p.add_argument("--protocol", type=Path, required=True)
    p.add_argument("--split", choices=tuple(COUNTS), required=True)
    p.add_argument("--index", type=int, required=True)
    p.add_argument("--map", type=Path, required=True)
    p.add_argument("--hapibd", type=Path, required=True)
    p.add_argument("--candidate", type=Path)
    p.add_argument("--baseline", type=Path)
    p.add_argument("--out", type=Path, required=True)
    for name in ("fit", "summarize"):
        p = sub.add_parser(name)
        p.add_argument("--protocol", type=Path, required=True)
        p.add_argument("--results", type=Path, required=True)
        p.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "protocol":
        result = make_protocol(args.out, args.map, args.hapibd, args.clustering_policy, args.selection_evidence)
    elif args.command == "run":
        result = run_cohort(args)
    elif args.command == "fit":
        result = fit(args)
    else:
        result = summarize(args)
    print(json.dumps({k: v for k, v in result.items() if k not in {"source_sha256", "bound_inputs", "inputs", "cohorts"}}, indent=2))
    if result.get("completed") is False:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
