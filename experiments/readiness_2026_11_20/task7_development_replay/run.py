"""Frozen development interoperability replay; no truth scoring or simulation."""
from __future__ import annotations

import hashlib
import importlib.metadata
import json
import os
import subprocess
import sys
import time
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SNAPSHOT = ROOT / 'source_snapshot'
sys.path[:0] = [str(SNAPSHOT / 'src'), str(SNAPSHOT)]

import numpy as np
import tskit
from ghist_single_sweep.export import write_biallelic_vcf
from ghist_single_sweep.scan import read_phased_biallelic_vcf
from ghist_submission.validate import validate_submission
from scripts.freeze_task7_sweep_bgs import verify_candidate


def digest(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            result.update(block)
    return result.hexdigest()


def atomic(path, value):
    temporary = path.with_name('.' + path.name + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def verify_inputs(protocol):
    for relative, expected in protocol['source_sha256'].items():
        if digest(SNAPSHOT / relative) != expected:
            raise ValueError('Frozen source changed: ' + relative)
    for name, expected in protocol['input_sha256'].items():
        if digest(ROOT / 'inputs' / name) != expected:
            raise ValueError('Frozen replay input changed: ' + name)
    for name, expected in protocol['runner_sha256'].items():
        if digest(ROOT / name) != expected:
            raise ValueError('Frozen replay runner changed: ' + name)
    if digest(protocol['input_tree']['remote_path']) != protocol['input_tree']['sha256']:
        raise ValueError('Retained tree differs from the fixed input hash')
    software = protocol['software']
    if digest(software['raisd_binary_path']) != software['raisd_binary_sha256']:
        raise ValueError('RAiSD binary changed')
    actual = {name: importlib.metadata.version(name) for name in software['versions']}
    if actual != software['versions'] or sys.version != software['python']:
        raise ValueError('Runtime versions differ from the frozen protocol')


def main():
    protocol = json.loads((ROOT / 'protocol.json').read_text())
    output = ROOT / 'output'
    if output.exists():
        raise FileExistsError('Preserve existing evidence; output already exists')
    verify_inputs(protocol)
    output.mkdir()
    record = {'state': 'RUNNING', 'protocol_sha256': digest(ROOT / 'protocol.json'),
              'slurm_job_id': os.environ.get('SLURM_JOB_ID'), 'planned_replays': 1,
              'scope': protocol['purpose'], 'submission_count': 0, 'release_ready': False,
              'scientifically_validated': False, 'seed': protocol['seed']}
    atomic(output / 'completion.json', record)
    started = time.monotonic()
    try:
        original = tskit.load(protocol['input_tree']['remote_path'])
        left = protocol['coordinate_transform']['original_left']
        right = protocol['coordinate_transform']['original_right']
        if original.num_samples != 50 or original.num_migrations != 0:
            raise ValueError('Original sample/migration structure differs from preflight')
        if (np.any(original.tables.edges.left < left)
                or np.any(original.tables.edges.right > right)
                or np.any(original.tables.sites.position < left)
                or np.any(original.tables.sites.position >= right)):
            raise ValueError('Original tree contains material outside the complete fixed segment')
        localized = original.keep_intervals([[left, right]], simplify=False).trim()
        if (localized.sequence_length != right - left
                or not np.array_equal(localized.samples(), original.samples())
                or not np.array_equal(localized.tables.sites.position,
                                      original.tables.sites.position - left)
                or not np.array_equal(localized.genotype_matrix(), original.genotype_matrix())):
            raise ValueError('Whole-segment translation changed genotypes, sites, samples or exposure')
        vcf = output / 'canonical_input.vcf'
        export = write_biallelic_vcf(localized, vcf, chromosome='7')
        reread = read_phased_biallelic_vcf(vcf)
        sample_nodes = list(localized.samples())
        groups = {}
        for node in sample_nodes:
            groups.setdefault(localized.node(node).individual, []).append(node)
        if any(individual >= 0 for individual in groups):
            if -1 in groups or any(len(nodes) != 2 for nodes in groups.values()):
                raise ValueError('Sample individuals are not diploid')
            sample_nodes = [node for nodes in groups.values() for node in nodes]
        expected_positions, expected_genotypes = [], []
        for variant in localized.variants(samples=sample_nodes):
            calls = variant.genotypes
            if (len(variant.alleles) == 2 and np.all((calls == 0) | (calls == 1))
                    and 0 < int(calls.sum()) < len(calls)):
                expected_positions.append(int(variant.site.position))
                expected_genotypes.append(calls)
        if (reread.contig_length != right - left
                or not np.array_equal(reread.positions, expected_positions)
                or not np.array_equal(reread.haplotypes, np.asarray(expected_genotypes))
                or export['accepted_biallelic_polymorphic_sites'] != len(reread.positions)):
            raise ValueError('Canonical tree export failed genotype/coordinate/length round trip')
        atomic(output / 'export.json', {**export, 'sample_node_order': [int(n) for n in sample_nodes],
                                        'tree_sha256': protocol['input_tree']['sha256'],
                                        'original_coordinate_origin': left,
                                        'vcf_sha256': digest(vcf), 'round_trip_verified': True})
        environment = dict(os.environ, PYTHONPATH=str(SNAPSHOT / 'src') + ':' + str(SNAPSHOT))
        command = [sys.executable, '-m', 'scripts.predict_task7_sweep_bgs',
                   '--vcf', str(vcf), '--accessible', str(ROOT / 'inputs/accessible.bed'),
                   '--dataset', 'testing', '--raisd', protocol['software']['raisd_binary_path'],
                   '--out-dir', str(output / 'prediction')]
        with (output / 'predictor.stdout').open('w') as stdout, (output / 'predictor.stderr').open('w') as stderr:
            subprocess.run(command, env=environment, cwd=ROOT, stdout=stdout, stderr=stderr, check=True)
        candidate = output / 'prediction/prediction.json'
        verified = verify_candidate(candidate)
        bed_validation = validate_submission(output / 'prediction/candidate.bed', 7, 'testing')
        expected_abstention = (verified['consensus_eligible'] is False
                               and verified['signal_ratio'] is None
                               and verified['distant_peak_score'] is None)
        verify_inputs(protocol)
        if digest(ROOT / 'protocol.json') != record['protocol_sha256']:
            raise ValueError('Frozen protocol changed during execution')
        if not expected_abstention:
            raise ValueError('Short segment unexpectedly acquired a distant-peak eligibility baseline')
        record.update(state='COMPLETED', canonical_round_trip=True, candidate_recomputed=True,
                      bed_validation=bed_validation, expected_short_segment_abstention=True,
                      callable_sites=verified['callable_sites'], callable_bases=verified['callable_bases'],
                      raisd_windows_total=verified['raisd_windows_total'],
                      raisd_windows_callable=verified['raisd_windows_callable'],
                      software=protocol['software'], no_localization_accuracy_scoring=True)
    except Exception as error:
        record.update(state='FAILED', failure=f'{type(error).__name__}: {error}')
        (output / 'failure.txt').write_text(traceback.format_exc())
        raise
    finally:
        record['elapsed_seconds'] = time.monotonic() - started
        record['output_sha256'] = {
            str(path.relative_to(output)): digest(path)
            for path in sorted(output.rglob('*')) if path.is_file() and path.name != 'completion.json'
        }
        atomic(output / 'completion.json', record)


if __name__ == '__main__':
    main()
