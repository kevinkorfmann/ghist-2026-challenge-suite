from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from ghist_growth_bgs.spectrum import fit_growth_spectrum
from ghist_inference.bootstrap import make_contract, replica_record, verify_inputs
from ghist_inference.checkpoint import atomic_json


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--blocks', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--replicate-start', type=int, required=True)
    parser.add_argument('--replicates', type=int, default=4)
    parser.add_argument('--planned-replicates', type=int, required=True,
                        help='Complete declared runner total, never the observed shard prefix')
    parser.add_argument('--projection', type=int, default=84)
    parser.add_argument('--restarts', type=int, default=12)
    parser.add_argument('--maxiter', type=int, default=400)
    parser.add_argument('--mutation-rate', type=float, default=3.36e-8)
    parser.add_argument('--seed', type=int, default=75194000)
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError('Preserve existing bootstrap artifact; use a fresh output')
    if args.replicates < 1 or args.replicate_start < 0 or args.replicate_start + args.replicates > args.planned_replicates:
        raise ValueError('Bootstrap task range must lie within its declared complete plan')
    settings = {'projection': args.projection, 'restarts': args.restarts, 'maxiter': args.maxiter,
                'mutation_rate': args.mutation_rate, 'convention': 'exponential', 'seed': args.seed,
                'resampling': 'within_contig_complete_blocks', 'fit_seed_stride': 1009}
    contract = make_contract(task=3, planned_replicates=args.planned_replicates,
                             blocks=args.blocks, settings=settings)
    with np.load(args.blocks, allow_pickle=False) as archive:
        spectra, exposure = archive['spectra'], archive['exposure_bp']
        contigs, starts, ends = archive['contig'], archive['start'], archive['end']
        block_size = int(archive['block_size'])
    verify_inputs(contract)
    usable = (ends - starts) == block_size
    groups = {contig: np.flatnonzero((contigs == contig) & usable) for contig in np.unique(contigs)}
    if any(len(indices) < 2 for indices in groups.values()):
        raise ValueError('each chromosome arm needs at least two full blocks')
    payload = {'bootstrap_contract': contract, 'state': 'RUNNING',
               'assigned_replicate_ids': list(range(args.replicate_start, args.replicate_start + args.replicates)),
               'projection': args.projection, 'block_size': block_size,
               'full_blocks_by_contig': {str(contig): len(indices) for contig, indices in groups.items()},
               'terminal_partial_blocks_excluded_from_resampling': True, 'results': []}
    atomic_json(args.out, payload)
    for replicate in payload['assigned_replicate_ids']:
        rng = np.random.default_rng(args.seed + replicate)
        selected = np.concatenate([rng.choice(indices, size=len(indices), replace=True) for indices in groups.values()])
        try:
            fit = fit_growth_spectrum(
                spectra[selected].sum(axis=0), sequence_length=float(exposure[selected].sum()),
                mutation_rate=args.mutation_rate, projection=args.projection, restarts=args.restarts,
                workers=1, maxiter=args.maxiter, seed=args.seed + 1009 * replicate)
            row = replica_record(replicate, selected.tolist(), fit, contract)
        except Exception as error:
            row = {'replicate': replicate, 'selected_blocks': selected.tolist(), 'state': 'FAILED',
                   'failure': f'{type(error).__name__}: {error}'}
        payload['results'].append(row)
        atomic_json(args.out, payload)
        print(f"replicate={replicate} state={row['state']}", flush=True)
    try:
        verify_inputs(contract)
        payload['state'] = 'COMPLETED' if all(row['state'] == 'COMPLETED' for row in payload['results']) else 'FAILED'
    except Exception as error:
        payload.update(state='FAILED', failure=f'{type(error).__name__}: {error}')
    atomic_json(args.out, payload)
    if payload['state'] != 'COMPLETED':
        raise RuntimeError('Bootstrap task failed; all attempted replicas retained and no intervals authorized')
    print(json.dumps({'state': payload['state'], 'replicates': len(payload['results'])}))


if __name__ == '__main__':
    main()
