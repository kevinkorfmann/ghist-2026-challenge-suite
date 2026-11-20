"""Independent, small complete-family operational check of frozen Task 10.

Generates new neutral coalescent founders and gene-drops them through one joint
18-person family, then runs the existing detected-IBD/branch/consensus workflow.
Oracle phase, a simple founder demography, and one small topology are explicit
limitations. This script never fits a model or accesses challenge inputs.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import itertools
import json
import os
import shutil
import subprocess
import sys
import time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from ghist_relatedness.provenance import bind_outputs


DEFINITION = "def complete 1 4\n2 0 3\n3 0 9\n4 2 9\n"
SOURCE_FILES = (
    "scripts/validate_task10_complete_pipeline.py",
    "scripts/prepare_hapibd_vcf.py",
    "scripts/summarize_hapibd.py",
    "scripts/summarize_hapibd_shape.py",
    "scripts/fit_human_family_relationships.py",
    "scripts/cluster_human_family_branches.py",
    "scripts/reconstruct_human_branch_pedigrees.py",
    "scripts/reconstruct_bonsai_pedigree.py",
    "scripts/consensus_human_family_pedigree.py",
    "src/ghist_relatedness/labels.py",
    "src/ghist_relatedness/provenance.py",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def save_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")


def save_result(output: Path, result: dict) -> None:
    save_json(output / "result.json", bind_outputs(
        output, result, protocol_name="preregistered_design.json"
    ))


def map_bounds(path: Path) -> dict[str, tuple[int, int, float]]:
    values: dict[str, list[tuple[int, float]]] = defaultdict(list)
    for line in path.read_text().splitlines():
        chrom, _, cm, bp = line.split()
        if chrom.isdigit() and 1 <= int(chrom) <= 22:
            values[chrom].append((int(bp), float(cm)))
    return {
        chrom: (rows[0][0], rows[-1][0], rows[-1][1] - rows[0][1])
        for chrom, rows in values.items()
    }


def founder_chromosome(job: tuple[str, int, int, float, int, str]) -> dict:
    """New 14-diploid-sample genealogy; no challenge or training genotypes."""
    import msprime
    import numpy as np

    chrom, start, end, cm, seed, filename = job
    length = end - start + 1
    ancestry = msprime.sim_ancestry(
        samples=14,
        population_size=10_000,
        sequence_length=length,
        recombination_rate=cm / (100 * length),
        random_seed=seed,
    )
    mutated = msprime.sim_mutations(
        ancestry, rate=1.25e-8, model=msprime.BinaryMutationModel(),
        random_seed=seed + 1,
    )
    retained = 0
    last = -1000
    with gzip.open(filename, "wt", compresslevel=1) as handle:
        for variant in mutated.variants():
            position = int(variant.site.position) + start
            counts = int(variant.genotypes.sum())
            # Fixed common-marker ascertainment and spacing, not selected on scores.
            if position - last < 500 or not 2 <= counts <= 26:
                continue
            genotypes = np.asarray(variant.genotypes).reshape(14, 2)
            calls = "\t".join(f"{a}|{b}" for a, b in genotypes)
            handle.write(f"{chrom}\t{position}\t.\tA\tC\t.\tPASS\t.\tGT\t{calls}\n")
            last = position
            retained += 1
    return {"chromosome": chrom, "variants": retained, "ancestry_seed": seed,
            "mutation_seed": seed + 1, "start_bp": start, "end_bp": end}


def generate_founders(bounds: dict, seed: int, root: Path, workers: int) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    jobs = [
        (chrom, *bounds[chrom], seed + 100 * int(chrom), str(root / f"chr{chrom}.vcf.gz"))
        for chrom in sorted(bounds, key=int)
    ]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        results = list(pool.map(founder_chromosome, jobs))
    path = root / "founders.vcf.gz"
    with gzip.open(path, "wt", compresslevel=1) as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        for chrom in sorted(bounds, key=int):
            handle.write(f"##contig=<ID={chrom},length={bounds[chrom][1] + 1}>\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
        handle.write("\t".join(f"founder_{i}" for i in range(14)) + "\n")
        for job in jobs:
            with gzip.open(job[-1], "rt") as source:
                shutil.copyfileobj(source, handle)
    save_json(root / "generation.json", {
        "founders": 14, "population_size": 10_000, "mutation_rate": 1.25e-8,
        "founder_recombination": "constant per chromosome using map mean",
        "marker_policy": "minor allele count >= 2 among 28 haplotypes; >=500 bp spacing",
        "chromosomes": results, "variants": sum(row["variants"] for row in results),
    })
    return path


def parents_from_fam(path: Path) -> dict[str, tuple[str, ...]]:
    parents = {}
    for line in path.read_text().splitlines():
        _, sample, father, mother, *_ = line.split()
        parents[sample] = tuple(parent for parent in (father, mother) if parent != "0")
    return parents


def ancestry_distances(sample: str, parents: dict) -> dict[str, int]:
    distances = {sample: 0}
    queue = [sample]
    while queue:
        current = queue.pop(0)
        for parent in parents[current]:
            distance = distances[current] + 1
            if parent not in distances or distances[parent] > distance:
                distances[parent] = distance
                queue.append(parent)
    return distances


def family_truth(parents: dict, identifiers: dict[str, str]) -> list[dict]:
    from ghist_relatedness.labels import accepted_label

    ancestors = {sample: ancestry_distances(sample, parents) for sample in identifiers}
    truth = []
    for left, right in itertools.combinations(sorted(identifiers), 2):
        common = set(ancestors[left]) & set(ancestors[right])
        if not common:
            raise ValueError("The chosen family must have a common ancestor for every pair")
        minimum = min(ancestors[left][a] + ancestors[right][a] for a in common)
        nearest = [a for a in common if ancestors[left][a] + ancestors[right][a] == minimum]
        distances = {(ancestors[left][a], ancestors[right][a]) for a in nearest}
        if len(distances) != 1:
            raise ValueError("This simple-tree truth oracle does not support unequal paths")
        up, down = next(iter(distances))
        label = accepted_label((up, down, len(nearest)), "human")
        if label is None:
            raise ValueError(f"Unsupported truth: {(up, down, len(nearest))}")
        a, b = sorted((identifiers[left], identifiers[right]))
        truth.append({"indID1": a, "indID2": b, "relation": label})
    return sorted(truth, key=lambda row: (row["indID1"], row["indID2"]))


def normalize_vcf(source: Path, destination: Path) -> tuple[dict, dict]:
    """Validate true phased calls and map synthetic IDs for Bonsai's numeric suffix."""
    identifiers = {}
    variants = 0
    heterozygotes = 0
    with gzip.open(source, "rt") as reader, gzip.open(destination, "wt", compresslevel=1) as writer:
        for line in reader:
            if line.startswith("#CHROM"):
                fields = line.rstrip().split("\t")
                identifiers = {sample: f"sample_{index}" for index, sample in enumerate(fields[9:])}
                fields[9:] = list(identifiers.values())
                writer.write("\t".join(fields) + "\n")
            elif line.startswith("#"):
                writer.write(line)
            else:
                fields = line.rstrip().split("\t")
                for value in fields[9:]:
                    genotype = value.split(":")[0]
                    if genotype not in {"0|0", "0|1", "1|0", "1|1"}:
                        raise ValueError(f"Oracle-phase assay rejects genotype {genotype!r}")
                    heterozygotes += genotype in {"0|1", "1|0"}
                variants += 1
                writer.write(line)
    if len(identifiers) != 18 or not variants or not heterozygotes:
        raise ValueError("Expected 18 samples and polymorphic phased genotype data")
    return identifiers, {"samples": len(identifiers), "variants": variants,
                         "phased_heterozygotes": heterozygotes, "unphased_heterozygotes": 0,
                         "missing_calls": 0, "phase_source": "Ped-sim --keep_phase (oracle phase)"}


def checked_rows(rows: list[dict], name: str) -> dict[tuple[str, str], dict]:
    output = {}
    for row in rows:
        pair = tuple(sorted((row["indID1"], row["indID2"])))
        if pair[0] == pair[1] or pair in output:
            raise ValueError(f"Self or duplicate pair in {name}: {pair}")
        output[pair] = row
    return output


def score_predictions(truth: list[dict], predictions: list[dict]) -> dict:
    """Missing predictions count as errors; no conditioning on observed IBD."""
    from ghist_relatedness.labels import HUMAN_LABELS

    meioses = {label: up + down for (up, down, _), label in HUMAN_LABELS.items()}
    expected = checked_rows(truth, "truth")
    actual = checked_rows(predictions, "predictions")
    if set(actual) - set(expected):
        raise ValueError("Predictions contain unknown pairs")
    if any(row["relation"] not in meioses for row in predictions):
        raise ValueError("Predictions contain an unsupported label")
    strata = {}
    correct = degree_correct = 0
    for pair, row in expected.items():
        label = row["relation"]
        prediction = actual.get(pair, {}).get("relation")
        match = prediction == label
        degree = prediction is not None and meioses[prediction] == meioses[label]
        correct += match
        degree_correct += degree
        counts = strata.setdefault(label, {"pairs": 0, "correct": 0, "meioses_correct": 0})
        counts["pairs"] += 1
        counts["correct"] += match
        counts["meioses_correct"] += degree
    return {"expected_pairs": len(expected), "predicted_pairs": len(actual),
            "missing_pairs": len(set(expected) - set(actual)),
            "all_pairs_present": set(expected) == set(actual),
            "exact_accuracy": correct / len(expected),
            "meioses_accuracy": degree_correct / len(expected), "strata": strata}


def frozen_pair_predictions(features: Path, candidate: Path, baseline: Path, output: Path) -> list[dict]:
    import joblib
    import numpy as np
    from fit_human_family_relationships import LABELS, matrix, read_observed

    def align_probabilities(source_classes, probabilities, target_classes):
        indices = {label: index for index, label in enumerate(source_classes)}
        if set(indices) != set(target_classes):
            raise ValueError("Frozen classifier label domain differs from the Task 10 contract")
        return probabilities[:, [indices[label] for label in target_classes]]

    records = read_observed(features)
    model = joblib.load(candidate)
    reference = joblib.load(baseline)
    classes = np.asarray(LABELS)
    weight = float(model["candidate_weight"])
    probabilities = weight * align_probabilities(
        model["candidate_model"].classes_,
        model["candidate_model"].predict_proba(matrix(records, model["candidate_feature_names"])), classes,
    ) + (1 - weight) * align_probabilities(
        reference["model"].classes_,
        reference["model"].predict_proba(matrix(records, reference["feature_names"])), classes,
    )
    if not np.isfinite(probabilities).all() or not np.allclose(probabilities.sum(axis=1), 1):
        raise ValueError("Invalid frozen-model probabilities")
    rows = [
        {"indID1": record["indID1"], "indID2": record["indID2"],
         "relation": str(classes[np.argmax(values)]),
         "probabilities": dict(zip(map(str, classes), map(float, values), strict=True))}
        for record, values in zip(records, probabilities, strict=True)
    ]
    save_json(output, rows)
    return rows


def run(args: argparse.Namespace) -> None:
    root, assets, output = args.project_root.resolve(), args.asset_root.resolve(), args.out_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    commands = []
    env = {**os.environ, "PYTHONPATH": f"{root / 'src'}:{root / 'scripts'}"}
    def script(name: str) -> list[str]:
        return [sys.executable, str(root / "scripts" / name)]

    def command(name: str, words: list) -> None:
        started = time.monotonic()
        with (output / f"{name}.log").open("w") as log:
            result = subprocess.run(list(map(str, words)), cwd=root, env=env, stdout=log, stderr=subprocess.STDOUT)
        commands.append({"stage": name, "command": list(map(str, words)),
                         "returncode": result.returncode, "seconds": time.monotonic() - started})
        save_json(output / "commands.json", commands)
        print(f"stage={name} returncode={result.returncode}", flush=True)
        if result.returncode:
            raise RuntimeError(f"Stage {name} failed; see {output / (name + '.log')}")

    map_path = assets / "external/ped-sim-map/refined_average.plink.map"
    simmap = assets / "external/ped-sim-map/refined_mf.simmap"
    candidate = assets / "artifacts/task10_large_family/clean_spectrum_v1/model.joblib"
    baseline = assets / "artifacts/task10_large_family/pair_model_min5/model.joblib"
    pedsim = assets / "external/ped-sim/ped-sim"
    interference = assets / "external/ped-sim/interfere/nu_p_campbell.tsv"
    hapibd = assets / "external/hap-ibd/hap-ibd.jar"
    # Save fixed design and provenance before generating or evaluating any data.
    provenance = {str(path): sha256(path) for path in
                  [map_path, simmap, candidate, baseline, pedsim, interference, hapibd]}
    provenance.update({str(root / path): sha256(root / path) for path in SOURCE_FILES})
    save_json(output / "preregistered_design.json", {
        "seed": args.seed, "pedsim_seed": args.seed + 9000, "definition": DEFINITION,
        "founders": 14, "genotyped_samples": 18, "expected_pairs": 153,
        "population_size": 10_000, "mutation_rate": 1.25e-8,
        "marker_policy": "minor allele count >=2; spacing >=500 bp",
        "phase": "oracle phase; no phasing algorithm executed",
        "hapibd": {"min_seed": 1, "min_extend": 0.2, "min_output": 2, "min_markers": 100, "max_gap": 1000},
        "branches": {"field": "relatedness_ge_20", "threshold": 0.12},
        "bonsai": {"threshold_cm": 10, "fast": [3, 1], "medium": [4, 2]},
        "consensus": "replicated topology and flat pair meiosis agreement; topology_boost=1",
        "success_criteria": "all stages complete, 153 unique valid-label pairs; accuracy descriptive only",
        "readiness_promotion": False, "model_fitting": False, "submission_count": 0,
        "hashes": provenance,
    })
    bounds = map_bounds(map_path)
    if len(bounds) != 22:
        raise ValueError("This operational assay requires all 22 autosomes")
    founders = generate_founders(bounds, args.seed, output / "founders", args.workers)
    definition = output / "family.def"
    definition.write_text(DEFINITION)
    command("ped_sim", [pedsim, "-d", definition, "-m", simmap, "-i", founders,
                        "-o", output / "family", "--intf", interference, "--seed", args.seed + 9000,
                        "--fam", "--keep_phase", "--err_rate", 0, "--miss_rate", 0])
    normalized = output / "family.normalized.vcf.gz"
    identifiers, vcf_audit = normalize_vcf(output / "family.vcf.gz", normalized)
    save_json(output / "identifiers.json", identifiers)
    save_json(output / "vcf_audit.json", vcf_audit)
    truth = family_truth(parents_from_fam(output / "family-everyone.fam"), identifiers)
    if len(truth) != 153:
        raise ValueError("Incomplete synthetic truth")
    save_json(output / "truth.json", truth)
    prepared = output / "prepared.vcf.gz"
    command("prepare", script("prepare_hapibd_vcf.py") + ["--input", normalized, "--output", prepared])
    prefix = output / "detected"
    command("hapibd", ["java", "-Xmx14g", "-jar", hapibd, f"gt={prepared}", f"map={map_path}",
                       f"out={prefix}", "min-seed=1.0", "min-extend=0.2", "min-output=2.0",
                       "min-markers=100", "max-gap=1000", f"nthreads={args.workers}"])
    ibd, hbd = output / "detected.ibd.gz", output / "detected.hbd.gz"
    command("summarize", script("summarize_hapibd.py") +
            ["--ibd", ibd, "--hbd", hbd, "--map", map_path, "--out-dir", output / "features"])
    features = output / "features/pair_shape_features.csv"
    command("shape", script("summarize_hapibd_shape.py") +
            ["--ibd", ibd, "--base", output / "features/pair_features.csv", "--out", features])
    predicted = frozen_pair_predictions(features, candidate, baseline, output / "pair_predictions.json")
    pair_score = score_predictions(truth, predicted)
    save_json(output / "pair_score.json", pair_score)
    # Fail closed on the production summarizer's omitted zero-IBD pairs. Never fill
    # them from truth or quietly report an accuracy conditional on detected pairs.
    if not pair_score["all_pairs_present"]:
        save_result(output, {"state": "FAILED_PAIR_COVERAGE", "pair_model": pair_score,
                                           "submission_ready": False, "submission_count": 0})
        raise RuntimeError("Production feature summarizer omitted one or more synthetic family pairs")
    cluster_file = output / "clusters.json"
    command("cluster", script("cluster_human_family_branches.py") +
            ["--pairs", features, "--threshold", 0.12, "--field", "relatedness_ge_20", "--out", cluster_file])
    for name, degree, types in [("fast", 3, 1), ("medium", 4, 2)]:
        command(name, script("reconstruct_human_branch_pedigrees.py") +
                ["--ibd", ibd, "--map", map_path, "--clusters", cluster_file, "--threshold-cm", 10,
                 "--max-add-degree", degree, "--max-rel-types", types, "--out", output / f"{name}.json"])
    command("consensus", script("consensus_human_family_pedigree.py") +
            ["--pair-predictions", output / "pair_predictions.json", "--bonsai", output / "fast.json",
             output / "medium.json", "--topology-boost", 1, "--require-flat-degree-agreement",
             "--out-dir", output / "consensus"])
    with (output / "consensus/candidate.csv").open() as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != ["indID1", "indID2", "relation"]:
            raise ValueError("Invalid candidate schema")
        consensus_score = score_predictions(truth, list(reader))
    if any(sha256(Path(path)) != digest for path, digest in provenance.items()):
        raise RuntimeError("Source, model or external tool changed during evaluation")
    result = {
        "state": "COMPLETED", "operational_contract_passed": consensus_score["all_pairs_present"],
        "independent_families": 1, "samples": 18, "truth_label_counts": dict(Counter(row["relation"] for row in truth)),
        "pair_model": pair_score, "consensus": consensus_score,
        "cluster_sizes": json.loads(cluster_file.read_text())["family_sizes"],
        "vcf_audit": vcf_audit, "submission_ready": False, "submission_count": 0,
        "limitations": ["oracle phase; no phasing inference", "one small 18-person family topology",
                        "only siblings, first and second cousins", "simple constant-Ne founder history",
                        "no genotype error or missingness", "not 229-person scaling or calibrated uncertainty"],
    }
    save_result(output, result)
    print(json.dumps(result, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--asset-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=2026112010)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        parser.error("Full-genome simulation must run as a scheduled Slurm job")
    run(args)


if __name__ == "__main__":
    main()
