from __future__ import annotations

import argparse
import hashlib
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
from pathlib import Path

import numpy as np

from ghist_bottleneck.sfs_baseline import (
    complete_block_spectra_from_vcf,
    fit_two_epoch_spectrum,
)
from ghist_inference.checkpoint import atomic_json, spectrum_digest


def validate_reference(reference, *, mutation_rate, sequence_length, folded):
    if bool(reference["folded"]) != folded:
        raise ValueError("reference fit and bootstrap folding settings disagree")
    if reference.get('sequence_length') != sequence_length:
        raise ValueError('reference fit and bootstrap sequence lengths disagree')
    if reference.get('callable_sequence_length', sequence_length) != sequence_length:
        raise ValueError('complete-block bootstrap requires full callable exposure')
    # Older fits did not persist mutation_rate, but theta/(4NL) identifies it.
    implied_rate = reference['theta'] / (4 * reference['ancestral_effective_size'] * sequence_length)
    if not np.isfinite(mutation_rate) or mutation_rate <= 0 or not np.isclose(implied_rate, mutation_rate, rtol=1e-10, atol=0):
        raise ValueError('reference fit and bootstrap mutation rates disagree')


def validate_input_reference(reference, vcf, blocks, metadata, *, allow_unverified=False):
    with Path(vcf).open('rb') as handle:
        vcf_hash = hashlib.file_digest(handle, 'sha256').hexdigest()
    spectrum_hash = spectrum_digest(blocks.sum(axis=0))
    reference_sfs = reference.get('sfs', {})
    if reference_sfs.get('projection') != metadata['projection']:
        raise ValueError('reference fit and bootstrap projections disagree')
    expected = {'input_vcf_sha256': reference.get('input_vcf_sha256'),
                'spectrum_sha256': reference_sfs.get('spectrum_sha256')}
    actual = {'input_vcf_sha256': vcf_hash, 'spectrum_sha256': spectrum_hash}
    for key, value in expected.items():
        if value is not None and value != actual[key]:
            raise ValueError(f'reference fit and bootstrap input disagree: {key}')
    missing = any(value is None for value in expected.values())
    if missing and not allow_unverified:
        raise ValueError('Reference lacks input provenance; regenerate the fit or explicitly use --allow-unverified-reference')
    return {**actual, 'reference_provenance': 'unverified_legacy_reference' if missing else 'verified_input_and_spectrum'}


def _fit_replica(arguments: tuple) -> dict:
    spectrum, mutation_rate, sequence_length, folded, initial, seed = arguments
    return fit_two_epoch_spectrum(
        spectrum,
        mutation_rate=mutation_rate,
        sequence_length=sequence_length,
        restarts=2,
        seed=seed,
        folded=folded,
        initial=tuple(initial),
    )


def record_replica(fit):
    """Retain estimates, but never turn failed optimizer termination into an interval."""
    result = {'fit': fit}
    if fit.get('optimizer', {}).get('converged') is not True:
        result['failure'] = 'Selected optimizer restart did not converge; intervals omitted'
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--fit", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--replicates", type=int, default=200)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--seed", type=int, default=83_771)
    parser.add_argument("--folded", action="store_true")
    parser.add_argument("--mutation-rate", type=float, default=3.42e-8)
    parser.add_argument("--sequence-length", type=int, default=100_000_000)
    parser.add_argument("--block-length", type=int, default=2_000_000)
    parser.add_argument('--allow-unverified-reference', action='store_true',
                        help='Explicitly permit a legacy fit without input hashes; result remains marked unverified')
    args = parser.parse_args()
    if args.replicates < 2 or args.workers < 1 or args.seed < 0:
        parser.error('at least two replicates, positive workers, and a nonnegative seed required')

    reference = json.loads(args.fit.read_text())
    validate_reference(reference, mutation_rate=args.mutation_rate,
                       sequence_length=args.sequence_length, folded=args.folded)
    blocks, metadata = complete_block_spectra_from_vcf(
        args.vcf,
        sequence_length=args.sequence_length,
        block_length=args.block_length,
        polarized=bool(reference["polarized"]),
    )
    if len(blocks) < 2:
        raise ValueError('bootstrap requires at least two independent resampling blocks')
    provenance = validate_input_reference(reference, args.vcf, blocks, metadata,
                                          allow_unverified=args.allow_unverified_reference)
    provenance['reference_fit_sha256'] = hashlib.sha256(args.fit.read_bytes()).hexdigest()
    rng = np.random.default_rng(args.seed)
    tasks = []
    initial = (reference["post_decline_fraction"], reference["scaled_time_2N"])
    for replica in range(args.replicates):
        selected = rng.integers(0, len(blocks), len(blocks))
        tasks.append(
            (
                blocks[selected].sum(axis=0),
                args.mutation_rate,
                args.sequence_length,
                args.folded,
                initial,
                args.seed + replica * 7_919,
            )
        )
    completed = {}
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        pending = {executor.submit(_fit_replica, task): index for index, task in enumerate(tasks)}
        for future in as_completed(pending):
            index = pending[future]
            try:
                completed[index] = record_replica(future.result())
            except Exception as error:
                completed[index] = {'failure': f'{type(error).__name__}: {error}'}
            atomic_json(args.out, dict(status='incomplete', requested_replicates=args.replicates,
                                      completed_replicates=completed, seed=args.seed, **provenance))
    failures = {index: result['failure'] for index, result in completed.items() if 'failure' in result}
    if failures:
        atomic_json(args.out, dict(status='failed', requested_replicates=args.replicates,
                                  completed_replicates=completed, seed=args.seed,
                                  **provenance,
                                  limitation='Intervals omitted: failed replicates cannot be silently dropped'))
        raise RuntimeError(f'{len(failures)} bootstrap fits failed; diagnostics saved to {args.out}')
    fits = [completed[index]['fit'] for index in range(args.replicates)]
    values = np.asarray(
        [[fit["generations"], fit["post_decline_fraction"], fit["ancestral_effective_size"]] for fit in fits]
    )
    bounded_fraction = np.clip(values[:, 1], 1e-6, 1 - 1e-6)
    latent = np.column_stack(
        (np.log(values[:, 0]), np.log(bounded_fraction / (1 - bounded_fraction)))
    )
    result = {
        **provenance,
        "status": "complete",
        "interval_kind": "block-bootstrap percentile; coverage has not been validated",
        "mutation_rate": args.mutation_rate,
        "sequence_length": args.sequence_length,
        "optimizer_diagnostics": [fit['optimizer'] for fit in fits],
        "converged_replicates": sum(fit['optimizer']['converged'] for fit in fits),
        "boundary_replicates": sum(bool(fit['optimizer_diagnostics']['near_boundary_parameter_indices']) for fit in fits),
        "replicates": args.replicates,
        "seed": args.seed,
        "folded": args.folded,
        "block_metadata": metadata,
        "columns": ["generations", "post_decline_fraction", "ancestral_effective_size"],
        "mean": values.mean(axis=0).tolist(),
        "median": np.median(values, axis=0).tolist(),
        "interval_95": np.quantile(values, [0.025, 0.975], axis=0).T.tolist(),
        "latent_covariance": np.cov(latent, rowvar=False).tolist(),
        "reference": {
            "generations": reference["generations"],
            "post_decline_fraction": reference["post_decline_fraction"],
            "ancestral_effective_size": reference["ancestral_effective_size"],
        },
        "values": values.tolist(),
    }
    atomic_json(args.out, result)
    print(json.dumps({key: result[key] for key in ("replicates", "mean", "median", "interval_95")}, indent=2))


if __name__ == "__main__":
    main()
