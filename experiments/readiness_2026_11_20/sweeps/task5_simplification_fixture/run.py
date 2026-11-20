"""Run the frozen pair and compare canonical biological/ancestry records."""
import hashlib
import json
from pathlib import Path
import subprocess
import time

import numpy as np
import pyslim
import tskit

HERE = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical(path):
    ts = tskit.load(path)
    individuals = [individual for individual in ts.individuals()
                   if individual.flags & pyslim.INDIVIDUAL_REMEMBERED]
    individuals.sort(key=lambda individual: individual.metadata['pedigree_id'])
    samples, identities = [], []
    for individual in individuals:
        nodes = sorted(individual.nodes, key=lambda node: ts.node(node).metadata['slim_id'])
        identities.append([individual.metadata['pedigree_id'], [ts.node(node).metadata['slim_id'] for node in nodes]])
        samples.extend(nodes)
    if len(individuals) != 12 or len(samples) != 24:
        raise ValueError('Expected exactly twelve remembered diploids')
    ts = ts.simplify(samples=samples, keep_input_roots=True, filter_populations=False)
    marginal = []
    for tree in ts.trees():
        clades = {node: tuple(sorted(tree.samples(node))) for node in tree.nodes()}
        rows = sorted((float(tree.time(node)), clades[node],
                       None if tree.parent(node) == tskit.NULL else
                       (float(tree.time(tree.parent(node))), clades[tree.parent(node)])) for node in tree.nodes())
        # Consecutive equal topologies can have implementation-specific breakpoints.
        if marginal and marginal[-1][2] == rows:
            marginal[-1][1] = float(tree.interval.right)
        else:
            marginal.append([float(tree.interval.left), float(tree.interval.right), rows])
    variants = []
    for variant in ts.variants():
        if len(variant.alleles) != 2 or np.any(variant.genotypes < 0):
            raise ValueError('Fixture expects complete biallelic selected-site genotypes')
        variants.append([float(variant.site.position), variant.genotypes.tolist()])
    return {'sample_identities': identities, 'marginal_ancestry': marginal, 'selected_genotypes': variants,
            'sequence_length': ts.sequence_length, 'sample_haplotypes': ts.num_samples}


def main():
    protocol_path = HERE / 'protocol.json'
    protocol_hash = digest(protocol_path)
    protocol = json.loads(protocol_path.read_text())
    locks = {HERE / name: value for name, value in protocol['source_sha256'].items()}
    locks.update({protocol_path: protocol_hash, Path(protocol['engine_path']): protocol['engine_sha256'],
                  Path(protocol['slim_path']): protocol['slim_sha256']})

    def verify():
        for path, expected in locks.items():
            if digest(path) != expected:
                raise ValueError(f'Frozen source/input drift: {path}')

    verify()
    output = HERE / 'output'
    output.mkdir(exist_ok=False)
    report = {'schema': 'paired_slim_simplification_software_result_v1', 'protocol_sha256': protocol_hash,
              'state': 'RUNNING', 'arms': {}, 'scientifically_validated': False,
              'scope': protocol['scope']}
    started = time.monotonic()
    representations = {}
    trajectories = {}
    try:
        for arm in ['default', 'interval20']:
            folder = output / arm
            folder.mkdir()
            begin = time.monotonic()
            with (folder / 'slim.log').open('w') as log:
                result = subprocess.run([protocol['slim_path'], '-s', str(protocol['seed']), str(HERE / f'{arm}.slim')],
                                        cwd=folder, stdout=log, stderr=subprocess.STDOUT,
                                        timeout=protocol['resource_limits']['per_arm_timeout_seconds'])
            entry = {'exit_code': result.returncode, 'elapsed_seconds': time.monotonic() - begin,
                     'log_sha256': digest(folder / 'slim.log')}
            report['arms'][arm] = entry
            if result.returncode != 0:
                raise RuntimeError(f'{arm} exited {result.returncode}; original log preserved')
            lines = (folder / 'slim.log').read_text().splitlines()
            trajectories[arm] = [line for line in lines if line.startswith('FIXTURE_TRAJECTORY')]
            if not trajectories[arm] or lines.count('FIXTURE_FORCED_RESTORE') != 1:
                raise ValueError(f'{arm} lacks trajectory or exactly one forced restore')
            representations[arm] = canonical(folder / 'result.trees')
            (folder / 'canonical.json').write_text(json.dumps(representations[arm], sort_keys=True) + '\n')
            entry.update(tree_sha256=digest(folder / 'result.trees'), canonical_sha256=digest(folder / 'canonical.json'),
                         trajectory_records=len(trajectories[arm]),
                         restore_debug_records=sum('restore()' in line for line in lines), forced_restore_records=1)
            verify()
        checks = {'trajectory': trajectories['default'] == trajectories['interval20']}
        checks.update({key: representations['default'][key] == representations['interval20'][key]
                       for key in representations['default']})
        report['comparisons'] = checks
        if not all(checks.values()):
            raise ValueError('Paired fixture comparison differs; mismatch preserved')
        report['state'] = 'PASSED'
    except Exception as error:
        report.update(state='FAILED', error=f'{type(error).__name__}: {error}')
        raise
    finally:
        report['elapsed_seconds'] = time.monotonic() - started
        (output / 'verification.json').write_text(json.dumps(report, indent=2) + '\n')
        print(json.dumps(report))


if __name__ == '__main__':
    main()
