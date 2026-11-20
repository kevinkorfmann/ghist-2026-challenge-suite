"""Prespecified target-length recovery with durable spectra and fit checkpoints.

Run once without --case to lock the protocol, then each of the 24 case indices.
Truth is used only for simulation and reporting, never optimizer initialization.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import tempfile
import time

import msprime
import numpy as np
import pysam
import tskit

from ghist_bottleneck.sfs_baseline import (
    complete_block_spectra_from_vcf, fit_two_epoch_spectrum, projected_sfs_from_vcf,
)
from ghist_bottleneck.dadi_reference import fit_two_epoch_dadi_spectrum
from ghist_bottleneck.optimization import optimizer_contract, validate_optimizer_contract
from ghist_inference.checkpoint import atomic_json, load_spectrum_checkpoint, save_spectrum_checkpoint
from scripts.validate_bottleneck_science import MU, RECOMBINATION, POINTS, panel_work


def protocol(length, seed, restarts):
    sources = [Path(__file__), Path('scripts/validate_bottleneck_science.py'),
               Path('src/ghist_bottleneck/sfs_baseline.py'),
               Path('src/ghist_bottleneck/dadi_reference.py'),
               Path('src/ghist_bottleneck/optimization.py'),
               Path('src/ghist_inference/validation.py'),
               Path('src/ghist_inference/checkpoint.py'),
               Path('src/ghist_inference/vcf.py')]
    return dict(
        schema='bottleneck_target_length_v3',
        optimizer_contracts={engine: optimizer_contract(engine) for engine in ('moments', 'dadi')}, sequence_length=length,
        block_length=length // 50, seed=seed, restarts=restarts,
        parameter_points=POINTS, replicates=2, cases=24, diploid_individuals=57,
        coalescent_ploidy=2, mutation_rate=MU, recombination_rate=RECOMBINATION,
        scenarios=['matched', 'ancestral_epoch'], folded=True,
        mutation_model='binary infinite sites; continuous positions',
        seed_allocation='disjoint consecutive ancestry/mutation/optimizer triples',
        fit_initialization='random log-uniform; no truth initialization',
        methods='moments on all cases; independent dadi diffusion on matched cases',
        vcf_checks='cases 0 and 12; binary alleles to A/C, floor(site)+1, exact contig length; compare all spectra',
        calibration_or_selection=False, public_answers_used=False,
        exposure=f'complete {length} bp simulation; all bases callable; no missing genotypes',
        limitations=[
            'Six diagnostic parameter points, not a random prior validation sample',
            'Two chromosomes per point/history cannot establish calibrated coverage',
            'Both numerical solvers assume the same two-epoch biological model',
            'No interval coverage claim; failures and incomplete cases remain in denominator',
            'Continuous sites may share an integer VCF position after rounding',
        ],
        source_sha256={str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        versions={p: importlib.metadata.version(p) for p in ['msprime', 'tskit', 'moments-popgen', 'dadi', 'numpy', 'scipy']},
    )


def vcf_contract(ts, length, expected_blocks, directory, identity, spectrum_sha256):
    tables = ts.dump_tables()
    tables.sites.packset_ancestral_state(['A'] * len(tables.sites))
    tables.mutations.packset_derived_state(['C' if m.derived_state == '1' else 'A' for m in ts.mutations()])
    dna = tables.tree_sequence()
    with tempfile.TemporaryDirectory(prefix='vcf-contract-', dir=directory) as temporary:
        path = Path(temporary) / 'simulation.vcf'
        with path.open('w') as handle:
            # tskit also transforms the sequence-length sentinel. Clamp that
            # sentinel, while every site still maps from [0,L) to [1,L].
            dna.write_vcf(handle, position_transform=lambda x: np.minimum(np.floor(x).astype(int) + 1, length))
        with pysam.VariantFile(str(path)) as vcf:
            if list(vcf.header.contigs) != ['1'] or vcf.header.contigs['1'].length != length:
                raise ValueError('Exported VCF has the wrong contig length')
        spectrum, metadata = projected_sfs_from_vcf(path, polarized=False)
        blocks, block_metadata = complete_block_spectra_from_vcf(
            path, sequence_length=length, block_length=length // 50, polarized=False,
        )
        np.testing.assert_array_equal(spectrum, expected_blocks.sum(axis=0))
        np.testing.assert_array_equal(blocks, expected_blocks)
    return dict(passed=True, identity=identity, spectrum_sha256=spectrum_sha256,
                contig_length=length, projected_sfs=metadata, block_sfs=block_metadata,
                diploid_individuals=dna.num_individuals, haplotypes=dna.num_samples)


def run_case(index, out, config, *, vcf_only=False):
    if not 0 <= index < 24:
        raise ValueError('case index must be in [0, 23]')
    point, replicate, scenario, length, seed, restarts = panel_work(
        2, config['sequence_length'], config['seed'], config['restarts'],
    )[index]
    directory = out / f'case_{index:02d}'
    directory.mkdir(parents=True, exist_ok=True)
    checkpoint = directory / 'spectrum.npz'
    identity = dict(case=index, point=point, replicate=replicate, scenario=scenario,
                    ancestry_seed=seed, mutation_seed=seed+1, optimizer_seed=seed+2,
                    protocol_sha256=hashlib.sha256((out / 'protocol.json').read_bytes()).hexdigest())
    if checkpoint.exists():
        blocks, metadata = load_spectrum_checkpoint(checkpoint, expected=identity)
    else:
        started = time.monotonic()
        generations, fraction, ne = POINTS[point]
        demography = msprime.Demography()
        demography.add_population(name='pop', initial_size=ne*fraction)
        demography.add_population_parameters_change(time=generations, initial_size=ne, population='pop')
        if scenario == 'ancestral_epoch':
            demography.add_population_parameters_change(time=4*generations, initial_size=1.8*ne, population='pop')
            demography.add_population_parameters_change(time=8*generations, initial_size=ne, population='pop')
        ancestry = msprime.sim_ancestry(
            samples={'pop': 57}, ploidy=2, demography=demography,
            sequence_length=length, recombination_rate=RECOMBINATION, random_seed=seed,
        )
        ts = msprime.sim_mutations(ancestry, rate=MU, model=msprime.BinaryMutationModel(),
                                  discrete_genome=False, random_seed=seed+1)
        blocks = ts.allele_frequency_spectrum(
            mode='site', polarised=True, span_normalise=False,
            windows=np.linspace(0, length, 51),
        )
        metadata = dict(**identity, truth=dict(generations=generations, post_decline_fraction=fraction,
                                              ancestral_effective_size=ne),
                        sites=ts.num_sites, trees=ts.num_trees,
                        simulation_seconds=time.monotonic()-started)
        # Persist before fitting, so a slow or interrupted optimizer cannot lose a simulation.
        if index in (0, 12):
            ts.dump(directory / 'vcf_pending.trees')
        metadata = save_spectrum_checkpoint(checkpoint, blocks, metadata)
        del ts, ancestry
    if index in (0, 12) and not (directory / 'vcf_contract.json').exists():
        ts = tskit.load(directory / 'vcf_pending.trees')
        check = vcf_contract(ts, length, blocks, directory, identity, metadata['spectrum_sha256'])
        atomic_json(directory / 'vcf_contract.json', check)
        (directory / 'vcf_pending.trees').unlink()
        del ts
    if index in (0, 12):
        check = json.loads((directory / 'vcf_contract.json').read_text())
        if (check.get('identity') != identity or check.get('spectrum_sha256') != metadata['spectrum_sha256']
                or check.get('passed') is not True or check.get('contig_length') != length):
            raise ValueError('VCF contract does not match this case, spectrum and protocol')
    if vcf_only:
        return
    spectrum = blocks.sum(axis=0)
    methods = {'moments': fit_two_epoch_spectrum}
    if scenario == 'matched':
        methods['dadi'] = fit_two_epoch_dadi_spectrum
    for name, fit in methods.items():
        target = directory / f'{name}.json'
        if target.exists():
            old = json.loads(target.read_text())
            if old.get('identity') != identity or old.get('spectrum_sha256') != metadata['spectrum_sha256']:
                raise ValueError('Fit checkpoint configuration mismatch')
            if 'fit' in old:
                validate_optimizer_contract(old['fit'], engine=name, expected=config['optimizer_contracts'][name], require_converged=False)
            continue
        started = time.monotonic()
        result = dict(identity=identity, truth=metadata['truth'], spectrum_sha256=metadata['spectrum_sha256'])
        try:
            estimate = fit(spectrum, mutation_rate=MU, sequence_length=length,
                           restarts=restarts, seed=seed+2, folded=True)
            validate_optimizer_contract(estimate, engine=name, expected=config["optimizer_contracts"][name], require_converged=False)
            result.update(fit=estimate, relative_error={k: estimate[k]/v-1 for k, v in metadata['truth'].items()})
        except Exception as error:
            result['failure'] = f'{type(error).__name__}: {error}'
        result['seconds'] = time.monotonic()-started
        atomic_json(target, result)
        print(json.dumps(dict(case=index, method=name, seconds=result['seconds'],
                              failed='failure' in result)), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--case', type=int)
    parser.add_argument('--sequence-length', type=int, default=100_000_000)
    parser.add_argument('--seed', type=int, default=2026112901)
    parser.add_argument('--restarts', type=int, default=8)
    parser.add_argument('--vcf-only', action='store_true', help='Revalidate export without repeating likelihood fits')
    args = parser.parse_args()
    if args.sequence_length < 50 or args.sequence_length % 50 or args.restarts < 1:
        parser.error('positive restarts and sequence length divisible by 50 required')
    panel_work(2, args.sequence_length, args.seed, args.restarts)
    config = protocol(args.sequence_length, args.seed, args.restarts)
    target = args.out / 'protocol.json'
    if target.exists():
        if json.loads(target.read_text()) != json.loads(json.dumps(config)):
            raise ValueError('Locked protocol mismatch; use a separate output directory')
    else:
        atomic_json(target, config)
    if args.case is not None:
        if args.vcf_only and args.case not in (0, 12):
            parser.error('--vcf-only requires case 0 or 12')
        run_case(args.case, args.out, config, vcf_only=args.vcf_only)


if __name__ == '__main__':
    main()
