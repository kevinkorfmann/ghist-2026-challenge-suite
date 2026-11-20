"""Experimental window-EM alternative; does not pass the realSFS release gates."""
from pathlib import Path
import re
import subprocess

import numpy as np

from .pipeline import audit_inputs, fit_spectrum, sha256, write_json


def read_winsfs(path: Path, sample_sizes: list[int]) -> np.ndarray:
    text = path.read_text()
    header = re.match(r'^#SHAPE=<([0-9/, ]+)>\s*\n', text)
    if not header:
        raise ValueError('missing winsfs shape header')
    shape = [int(s) for s in re.split(r'[,/ ]+', header[1])]
    if shape != [n+1 for n in sample_sizes]:
        raise ValueError('winsfs shape does not match samples')
    data = np.asarray(text[header.end():].split(), dtype=float)
    if data.size != np.prod(shape) or not np.all(np.isfinite(data)) or np.any(data < 0):
        raise ValueError('invalid winsfs values')
    return data


def make_safs(directory: Path, reference: Path, task: str, output: Path, angsd: Path,
              *, threads: int = 4, expected_length: int = 10_000_000):
    audit = audit_inputs(directory, reference, task, expected_length)
    output.mkdir(parents=True, exist_ok=False)
    write_json(output / 'inputs.json', audit)
    protocol = {'method': 'experimental_gl2_saf1_for_winsfs', 'source_sha256': sha256(Path(__file__)),
                'angsd_sha256': sha256(angsd), 'commands': []}
    indices = []
    for i, population in enumerate(audit['population_order']):
        listing = output / f'population_{i}.bams'
        listing.write_text(''.join(sample['path']+'\n' for sample in audit['samples'].values()
                                   if sample['population'] == population))
        prefix = output / f'population_{i}'
        command = [str(angsd), '-b', str(listing), '-ref', str(reference), '-anc', str(reference),
                   '-GL', '2', '-doSaf', '1', '-minMapQ', '30', '-minQ', '20', '-P', str(threads), '-out', str(prefix)]
        protocol['commands'].append(command)
        write_json(output / 'protocol.json', protocol)
        with (output / f'population_{i}.log').open('w') as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
        index = Path(str(prefix)+'.saf.idx')
        if not index.is_file() or index.stat().st_size == 0:
            raise RuntimeError('ANGSD produced no SAF index')
        indices.append(index)
    return audit, indices


def run_window_comparison(audit: dict, indices: list[Path], output: Path, winsfs: Path, *, threads: int = 4, max_epochs: int = 20):
    if max_epochs < 2 or threads < 1:
        raise ValueError('positive threads and at least two epochs required')
    output.mkdir(parents=True, exist_ok=False)
    protocol = {'method': 'experimental_winsfs0.7.0_published_defaults_vs_seed_continuation',
        'source_sha256': sha256(Path(__file__)), 'winsfs_sha256': sha256(winsfs),
        'saf_files': [{'path': str(member), 'sha256': sha256(member)} for p in indices
                      for member in [p, Path(str(p).removesuffix('.idx')+'.gz'), Path(str(p).removesuffix('.idx')+'.pos.gz')]],
        'sample_sizes': audit['sample_sizes'], 'commands': [],
        'tolerance': 1e-4, 'tolerance_units': 'winsfs sum of normalized block log likelihood changes; not equivalent to realSFS tolerance',
        'blocks': 500, 'window_size': 100, 'max_epochs': max_epochs,
        'release_gate_policy': 'v1 realSFS gates remain unpassed; separate development comparison only'}
    results, spectra = {}, {}
    for label, seed, initial in [('seed1', 761001, None), ('seed2', 761002, None),
                                  ('continued', 761001, output / 'seed1.unfolded.sfs')]:
        command = [str(winsfs), '-vv', '--seed', str(seed), '--threads', str(threads), '--blocks', '500',
                   '--window-size', '100', '--max-epochs', str(max_epochs), '--tolerance', '0.0001']
        if initial:
            command += ['--initial', str(initial)]
        command += [str(p) for p in indices]
        protocol['commands'].append(command)
        write_json(output / 'protocol.json', protocol)
        raw = output / f'{label}.unfolded.sfs'
        folded = output / f'{label}.folded.sfs'
        with raw.open('w') as out, (output / f'{label}.log').open('w') as log:
            subprocess.run(command, stdout=out, stderr=log, check=True)
        fold_command = [str(winsfs), 'view', '--fold', str(raw)]
        protocol['commands'].append(fold_command)
        write_json(output / 'protocol.json', protocol)
        with folded.open('w') as out:
            subprocess.run(fold_command, stdout=out, check=True)
        spectrum = read_winsfs(folded, audit['sample_sizes'])
        if not 0 < spectrum.sum() <= audit['length'] + 1:
            raise ValueError('invalid retained exposure from winsfs')
        estimate = fit_spectrum(spectrum, audit['sample_sizes'], audit['task'])
        estimate['callable_fraction'] = float(spectrum.sum()/audit['length'])
        estimate['winsfs_epochs'] = [int(n) for n in re.findall(r'Finished epoch ([0-9]+)', (output/f'{label}.log').read_text())]
        estimate['winsfs_early_stop_observed'] = bool(estimate['winsfs_epochs'] and estimate['winsfs_epochs'][-1] < max_epochs)
        write_json(output / f'{label}.estimate.json', estimate)
        results[label], spectra[label] = estimate, spectrum
    comparisons = {}
    for label in ['seed2', 'continued']:
        a, b = spectra['seed1'], spectra[label]
        comparisons[label] = {'frequency_L1': float(np.abs(a/a.sum()-b/b.sum()).sum()),
            'segregating_relative_L1': float(np.abs(a[1:]-b[1:]).sum()/a[1:].sum()),
            'target_relative_change': {k: results[label]['targets'][k]/v-1 for k, v in results['seed1']['targets'].items()}}
    result = {'task': audit['task'], 'estimates': {k: v['targets'] for k,v in results.items()},
        'comparisons': comparisons, 'optimizer_agreement': {k: v['optimizer_agreement'] for k,v in results.items()},
        'early_stop_observed': {k: v['winsfs_early_stop_observed'] for k,v in results.items()},
        'callable_fraction': {k: v['callable_fraction'] for k,v in results.items()},
        'scientifically_ready': False, 'release_gate_policy': protocol['release_gate_policy']}
    write_json(output / 'comparison.json', result)
    return result
