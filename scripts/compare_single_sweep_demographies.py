from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import stdpopsim

from ghist_single_sweep.scan import read_bed_intervals, read_phased_biallelic_vcf


CANDIDATES = [
    ("OutOfAfrica_3G09", "CEU"),
    ("OutOfAfrica_3G09", "CHB"),
    ("OutOfAfrica_2T12", "EUR"),
    ("OutOfAfricaExtendedNeandertalAdmixturePulse_3I21", "CEU"),
    ("OutOfAfrica_4J17", "CEU"),
    ("OutOfAfrica_4J17", "JPT"),
    ("Zigzag_1S14", "generic"),
    ("AmericanAdmixture_4B18", "ADMIX"),
    ("Africa_1T12", "AFR"),
]


def folded_summary(genotypes: np.ndarray, exposure: float) -> dict:
    sample_size = genotypes.shape[1]
    allele_counts = genotypes.sum(axis=1).astype(int)
    segregating = (allele_counts > 0) & (allele_counts < sample_size)
    allele_counts = allele_counts[segregating]
    minor = np.minimum(allele_counts, sample_size - allele_counts)
    spectrum = np.bincount(minor, minlength=sample_size // 2 + 1)[1:].astype(float)
    normalized = spectrum / max(spectrum.sum(), 1.0)
    pi = np.sum(2.0 * allele_counts * (sample_size - allele_counts)) / (
        sample_size * (sample_size - 1.0) * exposure
    )
    return {
        "sites": int(allele_counts.size),
        "pi": float(pi),
        "singleton_fraction": float(np.mean(minor == 1)),
        "intermediate_fraction": float(np.mean(minor >= sample_size * 0.2)),
        "folded_spectrum": spectrum.tolist(),
        "normalized_folded_spectrum": normalized.tolist(),
    }


def mean_summaries(summaries: list[dict]) -> dict:
    keys = ("sites", "pi", "singleton_fraction", "intermediate_fraction")
    result = {key: float(np.mean([item[key] for item in summaries])) for key in keys}
    result["normalized_folded_spectrum"] = np.mean(
        [item["normalized_folded_spectrum"] for item in summaries], axis=0
    ).tolist()
    return result


def distance(observed: dict, simulated: dict) -> dict:
    observed_sfs = np.asarray(observed["normalized_folded_spectrum"], dtype=float)
    simulated_sfs = np.asarray(simulated["normalized_folded_spectrum"], dtype=float)
    midpoint = (observed_sfs + simulated_sfs) / 2.0
    mask_observed = observed_sfs > 0
    mask_simulated = simulated_sfs > 0
    js = 0.5 * np.sum(observed_sfs[mask_observed] * np.log(observed_sfs[mask_observed] / midpoint[mask_observed]))
    js += 0.5 * np.sum(simulated_sfs[mask_simulated] * np.log(simulated_sfs[mask_simulated] / midpoint[mask_simulated]))
    components = {
        "sqrt_jensen_shannon": float(np.sqrt(max(js, 0.0))),
        "absolute_log_pi": float(abs(np.log(simulated["pi"] / observed["pi"]))),
        "absolute_singleton_fraction": float(
            abs(simulated["singleton_fraction"] - observed["singleton_fraction"])
        ),
        "absolute_intermediate_fraction": float(
            abs(simulated["intermediate_fraction"] - observed["intermediate_fraction"])
        ),
    }
    components["competition_distance"] = float(
        components["sqrt_jensen_shannon"]
        + 0.5 * components["absolute_log_pi"]
        + components["absolute_singleton_fraction"]
        + components["absolute_intermediate_fraction"]
    )
    return components


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--accessible-bed", type=Path, required=True)
    parser.add_argument("--exclude", type=int, nargs=2, default=[35_000_000, 41_000_000])
    parser.add_argument("--length", type=int, default=10_000_000)
    parser.add_argument("--replicates", type=int, default=3)
    parser.add_argument("--seed", type=int, default=51_000)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    variants = read_phased_biallelic_vcf(args.vcf)
    accessible = read_bed_intervals(args.accessible_bed, variants.contig)
    exposure = float(sum(right - left for left, right in accessible))
    exclude_mask = (variants.positions >= args.exclude[0]) & (variants.positions < args.exclude[1])
    exposure -= args.exclude[1] - args.exclude[0]
    observed = folded_summary(variants.haplotypes[~exclude_mask], exposure)

    species = stdpopsim.get_species("HomSap")
    contig = species.get_contig(
        length=args.length,
        mutation_rate=species.genome.mean_mutation_rate,
        recombination_rate=species.genome.mean_recombination_rate,
    )
    engine = stdpopsim.get_engine("msprime")
    results = []
    for candidate_index, (model_id, population) in enumerate(CANDIDATES):
        model = species.get_demographic_model(model_id)
        replicates = []
        for replicate in range(args.replicates):
            seed = args.seed + candidate_index * 100 + replicate
            ts = engine.simulate(
                model,
                contig,
                samples={population: 25},
                seed=seed,
            )
            matrix = ts.genotype_matrix().astype(np.uint8)
            replicates.append(folded_summary(matrix, ts.sequence_length))
        aggregate = mean_summaries(replicates)
        results.append(
            {
                "model": model_id,
                "population": population,
                "replicates": replicates,
                "aggregate": aggregate,
                "distance": distance(observed, aggregate),
            }
        )
    results.sort(key=lambda item: item["distance"]["competition_distance"])
    output = {
        "state": "COMPLETED",
        "observed_excluded_region": args.exclude,
        "observed": observed,
        "simulation_length": args.length,
        "replicates_per_candidate": args.replicates,
        "results": results,
        "selected": {
            "model": results[0]["model"],
            "population": results[0]["population"],
        },
        "submission_count": 0,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()

