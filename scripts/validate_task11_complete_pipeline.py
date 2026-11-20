"""Fixed, independent-family genotype smoke for Task 11; never a release gate.

The synthetic founder SNPs have no LD, and every observed family is a sibship.
This deliberately narrow test executes real detected-IBD/HBD stages; it cannot
establish matched dog-population or eight-label challenge accuracy.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
import gzip
import hashlib
from itertools import combinations
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import traceback

import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score

from ghist_relatedness.provenance import bind_outputs, verify_hashes

from scripts.cluster_task11_families import affinity_matrix, canonical_labels, cluster
from scripts.fit_task11_dog_relationships import feature_names, read_observed
from scripts.reconstruct_bonsai_pedigree import read_map


SEEDS = (2026112711, 2026112811)
MATING_TYPES = (
    "unrelated",
    "unrelated",
    "first_cousins",
    "first_cousins",
    "first_cousins",
    "siblings",
    "siblings",
    "siblings",
)
SIBSHIP_SIZE = 4
MARKER_SPACING_CM = 0.01
HBD_THRESHOLDS = (2.0, 5.0, 10.0)
ROOT = Path(__file__).resolve().parents[1]
SOURCE_FILES = (
    "scripts/validate_task11_complete_pipeline.py",
    "scripts/prepare_hapibd_vcf.py",
    "scripts/summarize_hapibd.py",
    "scripts/cluster_task11_families.py",
    "scripts/reconstruct_bonsai_pedigree.py",
    "scripts/fit_task11_dog_relationships.py",
    "scripts/reconstruct_task11_pedigrees.py",
    "scripts/consensus_task11_dog_pedigrees.py",
    "src/ghist_relatedness/labels.py",
    "src/ghist_relatedness/validation.py",
    "src/ghist_relatedness/provenance.py",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")


def gamete(parent: np.ndarray, cm: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Transmit founder ancestry using a Poisson process on the supplied dog map."""
    count = rng.poisson((cm[-1] - cm[0]) / 100)
    crossovers = np.sort(rng.uniform(cm[0], cm[-1], count))
    copy = (rng.integers(2) + np.searchsorted(crossovers, cm)) % 2
    return parent[copy, np.arange(cm.size)]


def family_ancestry(
    kind: str, cm: np.ndarray, rng: np.random.Generator, founder_offset: int
) -> tuple[np.ndarray, int]:
    next_haplotype = founder_offset

    def founder() -> np.ndarray:
        nonlocal next_haplotype
        result = np.repeat(
            np.arange(next_haplotype, next_haplotype + 2)[:, None], cm.size, axis=1
        )
        next_haplotype += 2
        return result

    def offspring(left: np.ndarray, right: np.ndarray) -> np.ndarray:
        return np.stack((gamete(left, cm, rng), gamete(right, cm, rng)))

    if kind == "unrelated":
        mother, father = founder(), founder()
    elif kind in {"siblings", "first_cousins"}:
        grandparent1, grandparent2 = founder(), founder()
        left = offspring(grandparent1, grandparent2)
        right = offspring(grandparent1, grandparent2)
        if kind == "siblings":
            mother, father = left, right
        else:
            mother = offspring(left, founder())
            father = offspring(right, founder())
    else:
        raise ValueError(f"Unknown mating type: {kind}")
    return np.stack(
        [offspring(mother, father) for _ in range(SIBSHIP_SIZE)]
    ), next_haplotype


def generate_cohort(map_path: Path, out: Path, seed: int) -> dict:
    """Generate all relatives jointly; no pair template or training data is read."""
    rng = np.random.default_rng(seed)
    maps = read_map(map_path)
    samples = [
        f"synthetic11_{index}" for index in range(len(MATING_TYPES) * SIBSHIP_SIZE)
    ]
    families = [
        samples[index * SIBSHIP_SIZE : (index + 1) * SIBSHIP_SIZE]
        for index in range(len(MATING_TYPES))
    ]
    hbd_cm = np.zeros(len(samples))
    markers = 0
    genome_cm = 0.0
    with gzip.open(out / "synthetic.vcf.gz", "wt") as destination:
        destination.write("##fileformat=VCFv4.2\n")
        destination.write(
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Known simulated phase">\n'
        )
        for chrom, (bp, _) in maps.items():
            destination.write(f"##contig=<ID={chrom},length={int(bp[-1])}>\n")
        destination.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
            + "\t".join(samples)
            + "\n"
        )
        for chrom, (map_bp, map_cm) in maps.items():
            requested_cm = np.linspace(
                map_cm[0],
                map_cm[-1],
                int(np.ceil((map_cm[-1] - map_cm[0]) / MARKER_SPACING_CM)) + 1,
            )
            bp = np.unique(np.interp(requested_cm, map_cm, map_bp).astype(int))
            cm = np.interp(bp, map_bp, map_cm)
            ancestries = []
            offset = 0
            for kind in MATING_TYPES:
                family, offset = family_ancestry(kind, cm, rng, offset)
                ancestries.append(family)
            ancestry = np.concatenate(ancestries)
            frequencies = rng.uniform(0.05, 0.5, size=cm.size)
            founder_alleles = rng.random((offset, cm.size)) < frequencies
            genotype = founder_alleles[ancestry, np.arange(cm.size)]
            is_hbd = ancestry[:, 0] == ancestry[:, 1]
            hbd_cm += (is_hbd[:, :-1] * np.diff(cm)).sum(axis=1)
            genome_cm += cm[-1] - cm[0]
            for index, position in enumerate(bp):
                calls = [
                    f"{int(pair[0])}|{int(pair[1])}" for pair in genotype[:, :, index]
                ]
                destination.write(
                    f"{chrom}\t{position}\t.\tA\tC\t.\tPASS\t.\tGT\t"
                    + "\t".join(calls)
                    + "\n"
                )
            markers += cm.size
    truth = {
        "seed": seed,
        "families": families,
        "samples": samples,
        "family_mating_types": list(MATING_TYPES),
        "relationship": "siblings",
        "true_within_family_pairs": len(families)
        * SIBSHIP_SIZE
        * (SIBSHIP_SIZE - 1)
        // 2,
        "genome_cm": float(genome_cm),
        "markers": markers,
        "marker_approximated_hbd_cm": dict(zip(samples, hbd_cm.tolist())),
        "expected_inbreeding": {
            "unrelated": 0.0,
            "first_cousins": 0.0625,
            "siblings": 0.25,
        },
        "independence": "All founder haplotypes, meioses, and whole families newly generated; no fitting.",
    }
    write_json(out / "truth.json", truth)
    return truth


def partition_metrics(truth: list[list[str]], predicted: list[list[str]]) -> dict:
    all_truth = [sample for family in truth for sample in family]
    all_predicted = [sample for family in predicted for sample in family]
    if len(set(all_predicted)) != len(all_predicted):
        raise ValueError("Duplicate sample in predicted family partition")
    if set(all_predicted) - set(all_truth):
        raise ValueError("Unknown sample in predicted family partition")
    truth_labels = {
        sample: index for index, family in enumerate(truth) for sample in family
    }
    labels = {
        sample: index for index, family in enumerate(predicted) for sample in family
    }
    true_pairs = {
        tuple(sorted(pair)) for family in truth for pair in combinations(family, 2)
    }
    predicted_pairs = {
        tuple(sorted(pair)) for family in predicted for pair in combinations(family, 2)
    }
    recovered = true_pairs & predicted_pairs
    exact = sum(set(family) in [set(item) for item in predicted] for family in truth)
    return {
        "truth_families": len(truth),
        "predicted_families": len(predicted),
        "truth_samples": len(all_truth),
        "covered_samples": len(all_predicted),
        "missing_samples": sorted(set(all_truth) - set(all_predicted)),
        "exactly_recovered_families": exact,
        "split_families": sum(
            len({labels[s] for s in family if s in labels}) > 1 for family in truth
        ),
        "merged_predicted_families": sum(
            len({truth_labels[s] for s in family}) > 1 for family in predicted
        ),
        "true_within_pairs": len(true_pairs),
        "predicted_within_pairs": len(predicted_pairs),
        "missed_within_pairs": len(true_pairs - predicted_pairs),
        "false_cross_family_pairs": len(predicted_pairs - true_pairs),
        "within_pair_recall": len(recovered) / len(true_pairs) if true_pairs else None,
        "within_pair_precision": len(recovered) / len(predicted_pairs)
        if predicted_pairs
        else None,
        "ari_with_missing_samples_as_distinct_singletons": float(
            adjusted_rand_score(
                [truth_labels[s] for s in all_truth],
                [
                    labels.get(s, len(predicted) + index)
                    for index, s in enumerate(all_truth)
                ],
            )
        ),
    }


def relationship_metrics(truth: dict, predictions: list[dict]) -> dict:
    predicted = {}
    for row in predictions:
        pair = tuple(sorted((row["indID1"], row["indID2"])))
        if pair in predicted:
            raise ValueError(f"Duplicate prediction: {pair}")
        if pair[0] == pair[1] or not set(pair).issubset(truth["samples"]):
            raise ValueError(f"Invalid prediction pair: {pair}")
        predicted[pair] = row["relation"]
    correct = covered = 0
    by_mating = defaultdict(lambda: {"correct": 0, "covered": 0, "pairs": 0})
    true_pairs = set()
    for family, kind in zip(
        truth["families"], truth["family_mating_types"], strict=True
    ):
        for pair in combinations(sorted(family), 2):
            true_pairs.add(pair)
            seen = pair in predicted
            right = predicted.get(pair) == truth["relationship"]
            correct += right
            covered += seen
            by_mating[kind]["pairs"] += 1
            by_mating[kind]["covered"] += seen
            by_mating[kind]["correct"] += right
    for values in by_mating.values():
        values["accuracy_counting_missing_as_wrong"] = (
            values["correct"] / values["pairs"]
        )
    return {
        "truth_pairs": len(true_pairs),
        "covered_true_pairs": covered,
        "missing_true_pairs": len(true_pairs) - covered,
        "correct": correct,
        "accuracy_counting_missing_as_wrong": correct / len(true_pairs),
        "extra_cross_family_predictions": len(set(predicted) - true_pairs),
        "by_parent_mating_type": dict(by_mating),
        "interpretation": "Sibling-only synthetic smoke; not balanced eight-label accuracy.",
    }


def run_command(name: str, command: list[str], out: Path, stages: list[dict]) -> None:
    with (
        (out / f"{name}.stdout.log").open("w") as stdout,
        (out / f"{name}.stderr.log").open("w") as stderr,
    ):
        result = subprocess.run(
            command, cwd=ROOT, stdout=stdout, stderr=stderr, check=False
        )
    stages.append({"stage": name, "command": command, "exit_code": result.returncode})
    if result.returncode:
        raise RuntimeError(
            f"{name} failed with exit {result.returncode}; see stage logs"
        )


def cluster_cohort(pair_path: Path, truth_sample_names: list[str], out: Path) -> dict:
    """Use genotype sample names only; no true partition is supplied to clustering."""
    pairs = pd.read_csv(pair_path)
    samples = sorted(truth_sample_names, key=lambda value: int(value.rsplit("_", 1)[1]))
    solutions = {}
    for threshold in (10, 20):
        affinity = affinity_matrix(pairs, samples, f"relatedness_ge_{threshold}")
        for seed in range(8):
            solutions[f"ge_{threshold}_seed_{seed}"] = canonical_labels(
                cluster(affinity, len(MATING_TYPES), seed), samples
            )
    reference = solutions["ge_10_seed_0"]
    minimum = min(
        float(adjusted_rand_score(reference, labels)) for labels in solutions.values()
    )
    families = [
        [
            sample
            for sample, label in zip(samples, reference, strict=True)
            if label == index
        ]
        for index in range(len(MATING_TYPES))
    ]
    result = {
        "state": "COMPLETED" if minimum == 1.0 else "UNSTABLE_DIAGNOSTIC_ONLY",
        "families": families,
        "minimum_adjusted_rand": minimum,
        "stability_runs": len(solutions),
        "family_count_supplied": len(MATING_TYPES),
        "reference_fixed_before_results": "ge_10_seed_0",
        "sample_roster_source": "synthetic VCF sample roster, including zero-IBD samples",
    }
    write_json(out / "clusters.json", result)
    return result


def run_pipeline(args: argparse.Namespace) -> dict:
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    paths = {
        "map": args.map.resolve(),
        "hapibd": args.hapibd.resolve(),
        "model": args.model.resolve(),
    }
    manifest = {
        "scope": "independent whole-family genotype operational smoke, not matched biological validation",
        "seed": args.seed,
        "fixed_seeds": list(SEEDS),
        "marker_spacing_cm": MARKER_SPACING_CM,
        "mating_types": list(MATING_TYPES),
        "sibship_size": SIBSHIP_SIZE,
        "hbd_thresholds": list(HBD_THRESHOLDS),
        "family_count_supplied": len(MATING_TYPES),
        "founder_model": "independent biallelic SNPs; common site frequency uniform(0.05,0.5)",
        "recombination": "Poisson crossovers on fixed dog average map",
        "phase": "known true simulation phase; missingness and genotype error zero",
        "fitting": "none; frozen pair forest and existing topology/consensus settings",
        "provenance": {
            key: {"path": str(path), "sha256": sha256(path)}
            for key, path in paths.items()
        },
        "sources": {name: sha256(ROOT / name) for name in SOURCE_FILES},
        "runtime": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "job_id": os.environ.get("SLURM_JOB_ID"),
        },
        "selection_policy": "No parameter adjustment from these outcomes; each entire cohort held out.",
        "failure_policy": "Record failing stage; missing true pairs remain in scoring denominators. An unstable clustering proceeds only as a fixed-reference diagnostic.",
        "submission_count": 0,
    }
    write_json(out / "protocol.json", manifest)
    bound_inputs = {
        **{str(path): manifest["provenance"][name]["sha256"] for name, path in paths.items()},
        **{str(ROOT / name): digest for name, digest in manifest["sources"].items()},
        str(out / "protocol.json"): sha256(out / "protocol.json"),
    }
    stages = []
    result = {
        "complete_task11_biological_validation": False,
        "submission_ready": False,
        "submission_count": 0,
        "stages": stages,
        "protocol_sha256": sha256(out / "protocol.json"),
    }
    truth = None
    try:
        truth = generate_cohort(paths["map"], out, args.seed)
        bound_inputs.update({
            str(out / name): sha256(out / name)
            for name in ("synthetic.vcf.gz", "truth.json")
        })
        stages.append(
            {
                "stage": "joint_pedigree_genotypes",
                "state": "COMPLETED",
                "vcf_sha256": sha256(out / "synthetic.vcf.gz"),
            }
        )

        def command(script: str, *values: object) -> list[str]:
            return [sys.executable, str(ROOT / "scripts" / script), *map(str, values)]

        run_command(
            "prepare",
            command(
                "prepare_hapibd_vcf.py",
                "--input",
                out / "synthetic.vcf.gz",
                "--output",
                out / "prepared.vcf.gz",
            ),
            out,
            stages,
        )
        run_command(
            "hapibd",
            [
                "java",
                "-Xmx10g",
                "-jar",
                str(paths["hapibd"]),
                f"gt={out / 'prepared.vcf.gz'}",
                f"map={paths['map']}",
                f"out={out / 'segments'}",
                "min-seed=1.0",
                "min-extend=0.2",
                "min-output=2.0",
                "min-markers=50",
                "max-gap=1000",
                "nthreads=4",
            ],
            out,
            stages,
        )
        run_command(
            "segment_summary",
            command(
                "summarize_hapibd.py",
                "--ibd",
                out / "segments.ibd.gz",
                "--hbd",
                out / "segments.hbd.gz",
                "--map",
                paths["map"],
                "--out-dir",
                out / "summary",
            ),
            out,
            stages,
        )
        clusters = cluster_cohort(
            out / "summary/pair_features.csv", truth["samples"], out
        )
        stages.append({"stage": "family_clustering", "state": clusters["state"]})
        result["family_partition"] = partition_metrics(
            truth["families"], clusters["families"]
        )
        result["cluster_stability"] = clusters["minimum_adjusted_rand"]
        hbd = pd.read_csv(out / "summary/hbd_features.csv").set_index("individual")
        result["hbd_by_parent_mating_type"] = {}
        for kind in sorted(set(MATING_TYPES)):
            members = [
                s
                for family, k in zip(truth["families"], MATING_TYPES, strict=True)
                if k == kind
                for s in family
            ]
            result["hbd_by_parent_mating_type"][kind] = {
                "samples": len(members),
                "expected_F": truth["expected_inbreeding"][kind],
                "mean_simulated_F_marker_approximation": float(
                    np.mean(
                        [
                            truth["marker_approximated_hbd_cm"][s] / truth["genome_cm"]
                            for s in members
                        ]
                    )
                ),
                "mean_detected_HBD_fraction_ge2cm": float(
                    np.mean(
                        [
                            hbd.loc[s, "inbreeding_fraction"] if s in hbd.index else 0.0
                            for s in members
                        ]
                    )
                ),
            }
        model = joblib.load(paths["model"])
        prediction_paths = []
        result["pair_model_by_hbd_threshold"] = {}
        for threshold in HBD_THRESHOLDS:
            features, pairs, family_indices = read_observed(
                out / "segments.ibd.gz",
                out / "segments.hbd.gz",
                paths["map"],
                out / "clusters.json",
                threshold,
            )
            if not np.isfinite(features).all():
                raise ValueError("Nonfinite HBD-masked features")
            probabilities = model.predict_proba(features)
            labels = model.classes_[np.argmax(probabilities, axis=1)]
            rows = [
                {
                    "indID1": pair[0],
                    "indID2": pair[1],
                    "family": family_indices[i],
                    "relation": str(labels[i]),
                    "probabilities": dict(
                        zip(model.classes_, probabilities[i].tolist())
                    ),
                    "features": dict(zip(feature_names(), features[i].tolist())),
                }
                for i, pair in enumerate(pairs)
            ]
            path = out / f"pair_predictions_hbd{int(threshold)}.json"
            write_json(path, rows)
            prediction_paths.append(path)
            result["pair_model_by_hbd_threshold"][str(threshold)] = (
                relationship_metrics(truth, rows)
            )
        stages.append({"stage": "hbd_mask_and_frozen_forest", "state": "COMPLETED"})
        run_command(
            "bonsai",
            command(
                "reconstruct_task11_pedigrees.py",
                "--ibd",
                out / "segments.ibd.gz",
                "--map",
                paths["map"],
                "--clusters",
                out / "clusters.json",
                "--genome-cm",
                truth["genome_cm"],
                "--out",
                out / "bonsai.json",
            ),
            out,
            stages,
        )
        run_command(
            "consensus",
            command(
                "consensus_task11_dog_pedigrees.py",
                "--pair-predictions",
                *prediction_paths,
                "--bonsai",
                out / "bonsai.json",
                "--out-dir",
                out / "consensus",
            ),
            out,
            stages,
        )
        result["consensus"] = relationship_metrics(
            truth, json.loads((out / "consensus/predictions.json").read_text())
        )
        verify_hashes(bound_inputs, context="Task 11 evaluation inputs changed")
        result["input_provenance_verified"] = True
        result["operational_stages_completed"] = True
        result["strict_cluster_stability_passed"] = (
            clusters["minimum_adjusted_rand"] == 1.0
        )
    except Exception as error:
        result["operational_stages_completed"] = False
        result["failure"] = {
            "type": type(error).__name__,
            "message": str(error),
            "traceback": traceback.format_exc(),
        }
        if truth is not None and "consensus" not in result:
            result["consensus"] = relationship_metrics(truth, [])
    result["remaining_limits"] = [
        "No independent dog founder haplotypes with matched LD/allele frequencies were available.",
        "Only sibling labels and four observed offspring per family; no challenge-sized heterogeneous topology.",
        "Known true phase and zero genotype error/missingness; no phasing validation.",
        "Eight families supplied to clustering; family-count inference is not evaluated.",
        "No posterior calibration, independent matched-population confirmation, or release certification.",
    ]
    result = bind_outputs(out, result, protocol_name="protocol.json")
    write_json(out / "result.json", result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--map", type=Path, required=True)
    parser.add_argument("--hapibd", type=Path, required=True)
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--seed", type=int, choices=SEEDS, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = run_pipeline(args)
    print(json.dumps(result, indent=2, allow_nan=False))
    if not result["operational_stages_completed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
