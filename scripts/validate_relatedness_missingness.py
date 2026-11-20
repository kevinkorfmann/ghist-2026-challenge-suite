"""Prespecified missing-call sensitivity on synthetic Task 11 cohorts.

These are repeated measurements of existing simulated families, not new
biological confirmation datasets. Neither policy is selected from the scores.
"""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path
import sys

import joblib
import numpy as np

from ghist_relatedness.provenance import bind_outputs, verify_hashes

from scripts.fit_task11_dog_relationships import read_observed
from scripts.prepare_hapibd_vcf import prepare_vcf
from scripts.validate_task11_complete_pipeline import (
    ROOT,
    cluster_cohort,
    partition_metrics,
    relationship_metrics,
    run_command,
    sha256,
    write_json,
)


RATES = (0.005, 0.02)
POLICIES = ("drop-site", "major")


def mask_genotypes(source: Path, target: Path, rate: float, seed: int) -> dict:
    if not 0 <= rate < 1:
        raise ValueError("Missingness rate must be between 0 and 1")
    rng = np.random.default_rng(seed)
    samples = records = masked_calls = incomplete_sites = 0
    with gzip.open(source, "rt") as src, gzip.open(target, "wt") as dst:
        for line in src:
            if line.startswith("#CHROM\t"):
                samples = len(line.rstrip().split("\t")) - 9
            if line.startswith("#"):
                dst.write(line)
                continue
            fields = line.rstrip("\r\n").split("\t")
            if samples < 1 or len(fields) != 9 + samples or fields[8] != "GT":
                raise ValueError("Expected the known synthetic GT-only cohort contract")
            if any(gt not in {"0|0", "0|1", "1|0", "1|1"} for gt in fields[9:]):
                raise ValueError(
                    "Baseline must have complete, phased biallelic genotypes"
                )
            missing = rng.random(samples) < rate
            records += 1
            masked_calls += int(missing.sum())
            incomplete_sites += int(missing.any())
            for index in np.flatnonzero(missing):
                fields[int(index) + 9] = ".|."
            dst.write("\t".join(fields) + "\n")
    if not records:
        raise ValueError("No synthetic variants")
    return {
        "seed": seed,
        "per_genotype_missing_probability": rate,
        "samples": samples,
        "records": records,
        "masked_genotypes": masked_calls,
        "sites_with_missing_genotypes": incomplete_sites,
        "expected_complete_site_fraction": (1 - rate) ** samples,
        "observed_complete_site_fraction": 1 - incomplete_sites / records,
        "process": "Independent whole-genotype masking; underlying alleles/phase unchanged",
    }


def run(args: argparse.Namespace) -> dict:
    baseline = args.baseline.resolve()
    truth = json.loads((baseline / "truth.json").read_text())
    original_protocol = json.loads((baseline / "protocol.json").read_text())
    original_result = json.loads((baseline / "result.json").read_text())
    baseline_hashes = original_result.get("output_hashes", {})
    # Scoring truth and the protocol must belong to the same recorded VCF run.
    # The cohort seed alone does not bind a pedigree, labels, or model identity.
    baseline_bindings = {
        str(baseline / name): baseline_hashes.get(name)
        for name in ("truth.json", "protocol.json", "synthetic.vcf.gz")
    }
    verify_hashes(baseline_bindings, context="Baseline cohort provenance")
    verify_hashes(
        {str(baseline / "protocol.json"): original_result.get("protocol_sha256")},
        context="Baseline result/protocol linkage",
    )
    if original_protocol.get("fixed_seeds") != [2026112711, 2026112811]:
        raise ValueError(
            "Requires the prespecified synthetic Task 11 diagnostic cohort"
        )
    if truth["seed"] != original_protocol["seed"]:
        raise ValueError("Synthetic cohort identity mismatch")
    source = baseline / "synthetic.vcf.gz"
    if sha256(source) != original_result["output_hashes"]["synthetic.vcf.gz"]:
        raise ValueError("Synthetic cohort changed after the baseline run")
    paths = {
        name: Path(original_protocol["provenance"][name]["path"])
        for name in ("map", "hapibd", "model")
    }
    for name, path in paths.items():
        if sha256(path) != original_protocol["provenance"][name]["sha256"]:
            raise ValueError(f"Baseline {name} artifact changed")
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    masking_seed = int(truth["seed"]) + 1_000_000 + RATES.index(args.rate)
    protocol = {
        "rates_fixed_before_outcomes": list(RATES),
        "policies_fixed_before_outcomes": list(POLICIES),
        "rate": args.rate,
        "masking_seed": masking_seed,
        "source_cohort_seed": truth["seed"],
        "source_vcf_sha256": sha256(source),
        "baseline_truth_sha256": sha256(baseline / "truth.json"),
        "baseline_protocol_sha256": sha256(baseline / "protocol.json"),
        "baseline_result_sha256": sha256(baseline / "result.json"),
        "fixed_hbd_threshold_cm": 5,
        "artifact_provenance": original_protocol["provenance"],
        "source_sha256": {
            name: sha256(ROOT / name)
            for name in (
                "scripts/validate_relatedness_missingness.py",
                "scripts/prepare_hapibd_vcf.py",
                "scripts/validate_task11_complete_pipeline.py",
                "scripts/summarize_hapibd.py",
                "scripts/fit_task11_dog_relationships.py",
                "scripts/cluster_task11_families.py",
                "src/ghist_relatedness/provenance.py",
            )
        },
        "scope": "Repeated missingness sensitivity of existing synthetic cohorts; no fitting or policy selection",
        "independent_new_biological_replicates": 0,
        "submission_ready": False,
    }
    write_json(out / "protocol.json", protocol)
    bound_inputs = {
        **baseline_bindings,
        str(baseline / "result.json"): protocol["baseline_result_sha256"],
        **{str(path): original_protocol["provenance"][name]["sha256"]
           for name, path in paths.items()},
        **{str(ROOT / name): digest for name, digest in protocol["source_sha256"].items()},
        str(out / "protocol.json"): sha256(out / "protocol.json"),
    }
    masked = out / "masked.vcf.gz"
    masking = mask_genotypes(source, masked, args.rate, masking_seed)
    bound_inputs[str(masked)] = sha256(masked)
    write_json(out / "masking.json", masking)
    results = {
        "protocol_sha256": sha256(out / "protocol.json"),
        "masking": masking,
        "source_cohort_seed": truth["seed"],
        "baseline_pair_model": original_result["pair_model_by_hbd_threshold"]["5.0"],
        "baseline_family_partition": original_result["family_partition"],
        "policies": {},
        "consensus": "Not reassessed: the separate baseline already fails consensus accuracy",
        "submission_ready": False,
    }
    model = joblib.load(paths["model"])
    for policy in POLICIES:
        policy_out = out / policy
        policy_out.mkdir()
        stages = []
        result = {"stages": stages, "status": "RUNNING"}
        results["policies"][policy] = result
        try:
            prepared = policy_out / "prepared.vcf.gz"
            result["preparation"] = prepare_vcf(masked, prepared, missing_policy=policy)
            write_json(policy_out / "preparation.json", result["preparation"])
            run_command(
                "hapibd",
                [
                    "java",
                    "-Xmx10g",
                    "-jar",
                    str(paths["hapibd"]),
                    f"gt={prepared}",
                    f"map={paths['map']}",
                    f"out={policy_out / 'segments'}",
                    "min-seed=1.0",
                    "min-extend=0.2",
                    "min-output=2.0",
                    "min-markers=50",
                    "max-gap=1000",
                    "nthreads=4",
                ],
                policy_out,
                stages,
            )
            run_command(
                "summary",
                [
                    sys.executable,
                    str(ROOT / "scripts/summarize_hapibd.py"),
                    "--ibd",
                    str(policy_out / "segments.ibd.gz"),
                    "--hbd",
                    str(policy_out / "segments.hbd.gz"),
                    "--map",
                    str(paths["map"]),
                    "--out-dir",
                    str(policy_out / "summary"),
                ],
                policy_out,
                stages,
            )
            clusters = cluster_cohort(
                policy_out / "summary/pair_features.csv", truth["samples"], policy_out
            )
            result["family_partition"] = partition_metrics(
                truth["families"], clusters["families"]
            )
            result["cluster_stability"] = clusters["minimum_adjusted_rand"]
            features, pairs, _ = read_observed(
                policy_out / "segments.ibd.gz",
                policy_out / "segments.hbd.gz",
                paths["map"],
                policy_out / "clusters.json",
                5.0,
            )
            if not np.isfinite(features).all():
                raise ValueError("Nonfinite detected-IBD features")
            labels = model.predict(features)
            rows = [
                {"indID1": pair[0], "indID2": pair[1], "relation": str(label)}
                for pair, label in zip(pairs, labels, strict=True)
            ]
            write_json(policy_out / "pair_predictions.json", rows)
            result["pair_model"] = relationship_metrics(truth, rows)
            verify_hashes(bound_inputs, context="Missingness evaluation inputs changed")
            result["input_provenance_verified"] = True
            result["status"] = "COMPLETED"
        except Exception as error:
            result["status"] = "FAILED"
            result["failure"] = {"type": type(error).__name__, "message": str(error)}
            result.setdefault("pair_model", relationship_metrics(truth, []))
        results = bind_outputs(out, results, protocol_name="protocol.json")
        write_json(out / "result.json", results)
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--rate", type=float, choices=RATES, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = run(args)
    print(json.dumps(result, indent=2, allow_nan=False))
    if any(row["status"] != "COMPLETED" for row in result["policies"].values()):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
