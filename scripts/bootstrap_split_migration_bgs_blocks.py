from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from ghist_split_migration_bgs.spectrum import fit_isolation_with_migration_spectrum
from ghist_inference.bootstrap import make_contract, replica_record, verify_inputs
from ghist_inference.checkpoint import atomic_json
from ghist_inference.freeze import require_optimizer_convergence


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--blocks', type=Path, required=True)
    parser.add_argument('--initial-fit', type=Path, required=True)
    parser.add_argument('--replicate', type=int, required=True)
    parser.add_argument('--planned-replicates', type=int, required=True,
                        help='Complete declared runner total, never the number of present files')
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--restarts', type=int, default=8)
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError('Preserve existing bootstrap artifact; use a fresh output')
    if not 0 <= args.replicate < args.planned_replicates:
        raise ValueError('Bootstrap replica must lie within its declared complete plan')
    settings = {'projection': [50, 40], 'restarts': args.restarts, 'maxiter': 600,
                'mutation_rate': 1.015e-8, 'resampling': 'within_contig_all_stored_blocks',
                'resampling_seed': 75198000, 'fit_seed': 75199000}
    contract = make_contract(task=4, planned_replicates=args.planned_replicates,
                             blocks=args.blocks, initial_fit=args.initial_fit, settings=settings)
    with np.load(args.blocks, allow_pickle=False) as archive:
        spectra, exposure, contigs = archive['spectra'], archive['exposure_bp'], archive['contig']
    initial = json.loads(args.initial_fit.read_text())['fit']
    verify_inputs(contract)
    require_optimizer_convergence(initial, 'bootstrap initial fit')
    if initial.get('projection') != [50, 40] or initial.get('mutation_rate') != settings['mutation_rate']:
        raise ValueError('Bootstrap initial fit projection/mutation rate differs from the protocol')
    scaled = initial['scaled_parameters']
    rng = np.random.default_rng(75198000 + args.replicate)
    selected = []
    for contig in np.unique(contigs):
        indices = np.flatnonzero(contigs == contig)
        selected.extend(rng.choice(indices, size=len(indices), replace=True).tolist())
    payload = {'bootstrap_contract': contract, 'state': 'RUNNING',
               'assigned_replicate_ids': [args.replicate], 'results': []}
    atomic_json(args.out, payload)
    try:
        fit = fit_isolation_with_migration_spectrum(
            np.sum(spectra[selected], axis=0), mutation_rate=settings['mutation_rate'],
            sequence_length=float(np.sum(exposure[selected])), projection=(50, 40),
            restarts=args.restarts, workers=1, maxiter=600, seed=75199000 + args.replicate,
            initial_scaled=(scaled['AFR_size_ratio'], scaled['EUR_size_ratio'],
                            scaled['split_duration_2N'], scaled['migration_2N']))
        row = replica_record(args.replicate, selected, fit, contract)
    except Exception as error:
        row = {'replicate': args.replicate, 'selected_blocks': selected, 'state': 'FAILED',
               'failure': f'{type(error).__name__}: {error}'}
    payload['results'].append(row)
    atomic_json(args.out, payload)
    try:
        verify_inputs(contract)
        payload['state'] = 'COMPLETED' if row['state'] == 'COMPLETED' else 'FAILED'
    except Exception as error:
        payload.update(state='FAILED', failure=f'{type(error).__name__}: {error}')
    atomic_json(args.out, payload)
    if payload['state'] != 'COMPLETED':
        raise RuntimeError('Bootstrap task failed; its attempted fit is retained and no intervals authorized')
    print(json.dumps({'state': payload['state'], 'replicate': args.replicate}))


if __name__ == '__main__':
    main()
