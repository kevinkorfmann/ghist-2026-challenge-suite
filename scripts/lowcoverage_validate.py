"""Independent synthetic BAM-to-demography validation, never challenge data.

Fixed round-valued designs were specified before any synthetic results. The
challenge testing-answer exposure is disclosed in the report; these designs
are not selected from those values. Two tasks use independent random seeds.
"""
from __future__ import annotations

import argparse
from array import array
import json
from pathlib import Path
import time

import msprime
import numpy as np
import pysam

from ghist_lowcoverage.pipeline import TASKS, MUTATION_RATE, run_pipeline, sha256, write_json


DESIGNS = {
    "bottleneck": {"ancestral_size": 10_000, "fraction": .2, "generations": 2_000, "seed": 731041},
    "split_isolation": {"ancestral_size": 10_000, "size1": 5_000, "size2": 15_000, "generations": 5_000, "seed": 731059},
}


def simulate_bams(task: str, output: Path, length: int, *, seed_offset: int = 0,
                  design_override: dict | None = None, depth_sigma: float = .25, quality: int = 30) -> dict:
    if task not in DESIGNS or length < 1000:
        raise ValueError("unknown task or sequence too short")
    design = dict(design_override or DESIGNS[task])
    if set(design) != set(DESIGNS[task]) or any(not np.isfinite(v) or v <= 0 for v in design.values()):
        raise ValueError("invalid explicit simulation design")
    if task == "bottleneck" and not 0 < design["fraction"] < 1:
        raise ValueError("simulation requires a population decline")
    if not np.isfinite(depth_sigma) or depth_sigma < 0 or quality not in (20, 30):
        raise ValueError("invalid simulated depth or quality profile")
    output.mkdir(parents=True, exist_ok=False)
    design["seed"] += seed_offset
    seed = design["seed"]
    model = msprime.Demography()
    if task == "bottleneck":
        model.add_population(name="pop1", initial_size=design["ancestral_size"] * design["fraction"])
        model.add_population_parameters_change(time=design["generations"], initial_size=design["ancestral_size"], population="pop1")
        targets = {"generations": design["generations"], "post_decline_fraction": design["fraction"]}
    else:
        model.add_population(name="pop1", initial_size=design["size1"])
        model.add_population(name="pop2", initial_size=design["size2"])
        model.add_population(name="ancestral", initial_size=design["ancestral_size"])
        model.add_population_split(time=design["generations"], derived=["pop1", "pop2"], ancestral="ancestral")
        targets = {"pop1_population_size": design["size1"], "pop2_population_size": design["size2"],
                   "generations_since_split": design["generations"]}
    error_probability = 10 ** (-quality / 10)
    record = {"task": task, "design": design, "truth": targets, "length": length, "mutation_rate": MUTATION_RATE,
              "recombination_rate": 1e-8, "read_length": 126, "fragment_length": 350,
              "error_probability": error_probability, "mapping": "perfect_known_alignment_no_BWA", "quality": quality,
              "individual_depth_lognormal_sigma": depth_sigma,
              "final_data_access": False, "testing_data_access": False,
              "limitation": "known alignment; no indels/mapping bias; finite-sites biallelic mutation model",
              "source_sha256": sha256(Path(__file__))}
    write_json(output / "design.json", record)
    ancestry = msprime.sim_ancestry(TASKS[task]["samples"], demography=model, sequence_length=length,
                                   recombination_rate=1e-8, ploidy=2, random_seed=seed)
    sequence = msprime.sim_mutations(ancestry, rate=MUTATION_RATE, random_seed=seed + 1, model=msprime.BinaryMutationModel())
    sequence.dump(str(output / "truth.trees"))
    rng = np.random.default_rng(seed + 2)
    bases = np.array(list("ACGT"), dtype="S1")
    reference_codes = rng.integers(0, 4, size=length, dtype=np.uint8)
    reference = bases[reference_codes].tobytes().decode()
    ref_path = output / "reference.fasta"
    with ref_path.open("w") as ref:
        ref.write(">1\n")
        for start in range(0, length, 80):
            ref.write(reference[start:start + 80] + "\n")
    # ANGSD compares integer-second mtimes with >=; same-second FAI creation fails.
    time.sleep(max(0.0, int(ref_path.stat().st_mtime) + 1.02 - time.time()))
    pysam.faidx(str(ref_path))
    variants = list(sequence.variants())
    positions = np.array([int(v.site.position) for v in variants])
    genotypes = np.array([v.genotypes for v in variants], dtype=np.uint8).reshape(len(variants), sequence.num_samples)
    sample_populations = [pop for pop, count in TASKS[task]["samples"].items() for _ in range(count)]
    sample_names = [f"synthetic_{i:02}" for i in range(len(sample_populations))]
    (output / "popfile.txt").write_text("".join(f"{sample} {pop}\n" for sample, pop in zip(sample_names, sample_populations)))
    multipliers = np.exp(rng.normal(0, depth_sigma, len(sample_names)))
    multipliers /= multipliers.mean()
    counts = []
    for index, sample in enumerate(sample_names):
        # Each sample has its own independently spawned random stream.
        srng = np.random.default_rng(np.random.SeedSequence([seed, index, 93829]))
        haps = np.repeat(reference_codes[None, :], 2, axis=0)
        for hap in range(2):
            active = genotypes[:, 2 * index + hap] == 1
            at = positions[active]
            haps[hap, at] = (reference_codes[at] + 1) % 4
        count = int(srng.poisson(TASKS[task]["depth"] * multipliers[index] * length / 252))
        header = {"HD": {"VN": "1.6", "SO": "unsorted"}, "SQ": [{"SN": "1", "LN": length}],
                  "RG": [{"ID": sample, "SM": sample}]}
        unsorted = output / f"{sample}.unsorted.bam"
        with pysam.AlignmentFile(str(unsorted), "wb", header=header) as bam:
            for pair in range(count):
                start = int(srng.integers(0, length - 350 + 1))
                hap = int(srng.integers(0, 2))
                for mate, pos in enumerate([start, start + 350 - 126]):
                    codes = haps[hap, pos:pos + 126].copy()
                    errors = srng.random(126) < error_probability
                    codes[errors] = (codes[errors] + srng.integers(1, 4, size=int(errors.sum()))) % 4
                    read = pysam.AlignedSegment()
                    read.query_name = f"{sample}_{pair}"
                    read.query_sequence = bases[codes].tobytes().decode()
                    read.query_qualities = array("B", [quality] * 126)
                    read.flag = 99 if mate == 0 else 147
                    read.reference_id = 0
                    read.reference_start = pos
                    read.mapping_quality = 60
                    read.cigar = [(0, 126)]
                    read.next_reference_id = 0
                    read.next_reference_start = start + 350 - 126 if mate == 0 else start
                    read.template_length = 350 if mate == 0 else -350
                    read.set_tag("RG", sample)
                    bam.write(read)
        sorted_bam = output / f"{sample}.bam"
        pysam.sort("-o", str(sorted_bam), str(unsorted))
        pysam.index(str(sorted_bam))
        unsorted.unlink()
        counts.append({"sample": sample, "pairs": count, "realized_mean_depth": count * 252 / length,
                       "sha256": sha256(sorted_bam)})
    spectrum = sequence.allele_frequency_spectrum(sample_sets=[sequence.samples(population=sequence.population(p).id)
                                      for p in range(len(TASKS[task]["samples"]))], mode="site", polarised=False, span_normalise=False)
    np.save(output / "truth_folded_sfs.npy", spectrum)
    record["samples"] = counts
    record["reference_sha256"] = sha256(ref_path)
    record["trees_sha256"] = sha256(output / "truth.trees")
    record["design_sha256"] = sha256(output / "design.json")
    write_json(output / "simulation.json", record)
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task", choices=DESIGNS, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--length", type=int, default=1_000_000)
    parser.add_argument("--seed-offset", type=int, default=0)
    parser.add_argument("--angsd", type=Path, required=True)
    parser.add_argument("--realsfs", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--restarts", type=int, default=8)
    parser.add_argument("--bwa", type=Path, help="Optional BWA-MEM realignment of generated read pairs")
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    args.output.mkdir(parents=True)
    start = time.time()
    truth = simulate_bams(args.task, args.output / "synthetic", args.length, seed_offset=args.seed_offset)
    input_directory = args.output / "synthetic"
    if args.bwa:
        from ghist_lowcoverage.remap import remap_synthetic_bams
        remap_synthetic_bams(input_directory, args.output / "remapped", args.bwa, args.threads)
        input_directory = args.output / "remapped"
    estimate = run_pipeline(input_directory, input_directory / "reference.fasta", args.task,
                            args.output / "inference", args.angsd, args.realsfs, threads=args.threads,
                            expected_length=args.length, restarts=args.restarts)
    relative = {key: estimate["targets"][key] / value - 1 for key, value in truth["truth"].items()}
    record = {"task": args.task, "length": args.length, "seed_offset": args.seed_offset, "truth": truth["truth"],
              "estimate": estimate["targets"], "relative_error": relative,
              "diagnostic_rrmse": float(np.sqrt(np.mean(np.square(list(relative.values()))))),
              "metric_note": "diagnostic only; official current scoring formula not independently verified",
              "elapsed_seconds": time.time() - start, "scientifically_ready": False,
              "source_sha256": sha256(Path(__file__)), "simulation_sha256": sha256(args.output / "synthetic/simulation.json"),
              "estimate_sha256": sha256(args.output / "inference/estimate.json"),
              "testing_data_access": False, "final_data_access": False,
              "mapping": "bwa_mem" if args.bwa else "known_correct_alignment",
              "numerically_valid": estimate["numerically_valid"],
              "remap_protocol_sha256": sha256(args.output / "remapped/remap_protocol.json") if args.bwa else None}
    write_json(args.output / "result.json", record)
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    main()
