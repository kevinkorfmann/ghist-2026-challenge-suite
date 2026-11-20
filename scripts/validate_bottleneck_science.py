"""Independent diploid recovery panel; no observed data or published answers."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import importlib.metadata
import json
from pathlib import Path

import msprime

from ghist_bottleneck.sfs_baseline import fit_two_epoch_spectrum
from ghist_bottleneck.dadi_reference import fit_two_epoch_dadi_spectrum


# Prespecified diagnostic points spanning recent/old and strong/weak declines.
# These are not a random draw from the challenge prior or a calibration set.
POINTS = [(600, .12, 8000), (1900, .35, 17000), (8500, .50, 32000),
          (16000, .72, 44000), (28000, .89, 21000), (4000, .94, 11000)]
MU = 3.42e-8
RECOMBINATION = 1.87e-8


def panel_work(replicates: int, length: int, seed: int, restarts: int) -> list[tuple]:
    """Allocate a disjoint ancestry/mutation/optimizer seed triple to each case."""
    count = len(POINTS) * replicates * 2
    if replicates < 1 or seed < 1 or seed + 3 * count - 1 >= 2**32:
        raise ValueError("Panel seeds must be positive distinct 32-bit integers")
    cases = [(point, rep, scenario) for point in range(len(POINTS))
             for rep in range(replicates) for scenario in ('matched', 'ancestral_epoch')]
    return [(*case, length, seed + 3 * index, restarts) for index, case in enumerate(cases)]


def validate_reused_cases(cases: list[dict], replicates: int) -> None:
    expected = {(point, rep, scenario) for point in range(len(POINTS))
                for rep in range(replicates) for scenario in ('matched', 'ancestral_epoch')}
    identities = [(case['point'], case['replicate'], case['scenario']) for case in cases]
    seeds = [case['seed'] + offset for case in cases for offset in range(3)]
    if len(identities) != len(expected) or set(identities) != expected:
        raise ValueError("Reused panel has duplicate or missing cases")
    if len(set(seeds)) != len(seeds) or min(seeds) < 1 or max(seeds) >= 2**32:
        raise ValueError("Reused panel has colliding or invalid simulation seeds")


def run_case(args: tuple) -> dict:
    point, replicate, scenario, length, case_seed, restarts = args
    generations, fraction, ne = POINTS[point]
    demography = msprime.Demography()
    demography.add_population(name="pop", initial_size=ne * fraction)
    demography.add_population_parameters_change(time=generations, initial_size=ne, population="pop")
    if scenario == "ancestral_epoch":
        demography.add_population_parameters_change(time=4*generations, initial_size=1.8*ne, population="pop")
        demography.add_population_parameters_change(time=8*generations, initial_size=ne, population="pop")
    ancestry = msprime.sim_ancestry(
        samples=[msprime.SampleSet(114, population="pop", ploidy=1)],
        ploidy=2, demography=demography, sequence_length=length,
        recombination_rate=RECOMBINATION, random_seed=case_seed,
    )
    ts = msprime.sim_mutations(ancestry, rate=MU, model=msprime.BinaryMutationModel(),
                              discrete_genome=False, random_seed=case_seed+1)
    spectrum = ts.allele_frequency_spectrum(mode="site", polarised=True, span_normalise=False)
    result = dict(point=point, replicate=replicate, scenario=scenario, seed=case_seed,
                  ancestry_seed=case_seed, mutation_seed=case_seed+1, optimizer_seed=case_seed+2,
                  truth=dict(generations=generations, post_decline_fraction=fraction, ancestral_effective_size=ne),
                  segregating_sites=ts.num_sites, spectrum=spectrum.tolist(), methods={})
    return fit_record((result, length, restarts))


def fit_record(args: tuple) -> dict:
    result, length, restarts = args
    result = {**result, "methods": {}}
    spectrum = result["spectrum"]
    methods = {"moments": fit_two_epoch_spectrum}
    if result["scenario"] == "matched":
        methods["dadi"] = fit_two_epoch_dadi_spectrum
    for name, fit in methods.items():
        try:
            estimate = fit(spectrum, mutation_rate=MU, sequence_length=length,
                           restarts=restarts, seed=result["seed"]+2, folded=True)
            result["methods"][name] = dict(
                fit=estimate,
                relative_error={key: estimate[key]/value-1 for key, value in result["truth"].items()},
            )
        except Exception as error:
            result["methods"][name] = {"failure": f"{type(error).__name__}: {error}"}
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--replicates", type=int, default=2)
    parser.add_argument("--sequence-length", type=int, default=10_000_000)
    parser.add_argument("--restarts", type=int, default=8)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--seed", type=int, default=202611201)
    parser.add_argument("--spectra-from", type=Path,
                        help="Refit a complete panel without resimulating chromosomes")
    args = parser.parse_args()
    if min(args.replicates, args.sequence_length, args.restarts, args.workers, args.seed) <= 0:
        parser.error("all numeric settings must be positive")
    work = panel_work(args.replicates, args.sequence_length, args.seed, args.restarts)
    protocol = dict(
        seed_allocation="consecutive disjoint ancestry/mutation/optimizer triples",
        coalescent_ploidy=2, haplotypes=114, parameter_points=POINTS,
        replicates=args.replicates, sequence_length=args.sequence_length,
        mutation_rate=MU, recombination_rate=RECOMBINATION, seed=args.seed,
        restarts=args.restarts, fit_initialization="random; no truth initialization",
        heldout_unit="independent chromosome; repeated parameter points are reported together",
        calibration_or_selection=False, public_answers_used=False,
        missingness="complete calls only; missing-call robustness not established",
        uncertainty_coverage="not assessed; bootstrap intervals are not calibration evidence",
        rate_sensitivity="At fixed SFS, inferred N and time scale inversely with assumed mutation rate; decline fraction is invariant. This structural sensitivity is not a recovery experiment.",
        limitation="10 Mb diagnostic panel is not full-prior or target 100 Mb validation; dadi and moments share the same two-epoch biological model",
        source_sha256={str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in [
            Path(__file__), Path('src/ghist_bottleneck/sfs_baseline.py'), Path('src/ghist_bottleneck/dadi_reference.py')]},
        versions={p:importlib.metadata.version(p) for p in ['msprime','tskit','moments-popgen','dadi','numpy']},
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    existing = None
    if args.spectra_from:
        existing = json.loads(args.spectra_from.read_text())
        for key in ('coalescent_ploidy', 'haplotypes', 'parameter_points', 'replicates',
                    'sequence_length', 'mutation_rate', 'recombination_rate', 'seed'):
            if json.dumps(existing['protocol'][key]) != json.dumps(protocol[key]):
                raise ValueError(f"Reused panel disagrees on {key}")
        if not existing['complete'] or len(existing['cases']) != len(POINTS)*args.replicates*2:
            raise ValueError("Cannot refit an incomplete panel")
        validate_reused_cases(existing['cases'], args.replicates)
        protocol['seed_allocation'] = existing['protocol'].get('seed_allocation', 'legacy arithmetic; uniqueness verified before reuse')
        protocol['reused_spectra_sha256'] = hashlib.sha256(args.spectra_from.read_bytes()).hexdigest()
    args.out.with_suffix('.protocol.json').write_text(json.dumps(protocol, indent=2)+'\n')
    records = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        results = pool.map(fit_record, [(case, args.sequence_length, args.restarts)
                                      for case in existing['cases']]) if existing else pool.map(run_case, work)
        for result in results:
            records.append(result)
            args.out.write_text(json.dumps(dict(protocol=protocol, cases=records, complete=len(records)==len(work)), indent=2)+'\n')
            print(f"{len(records)}/{len(work)} point={result['point']} scenario={result['scenario']}", flush=True)


if __name__ == '__main__':
    main()
