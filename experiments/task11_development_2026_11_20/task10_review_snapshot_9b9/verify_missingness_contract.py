"""Verify explicit missing-data rejection; never a robustness experiment."""
import argparse
import gzip
import json
import os
from contextlib import ExitStack
from pathlib import Path

import numpy as np

from scripts.develop_task10_detected_ibd import save_json, sha256
from scripts.prepare_hapibd_vcf import prepare_vcf


def mask_genotypes(source, destinations, seeds):
    """Nested MCAR genotype masks, invariant to input sample-column order."""
    generators = {chrom: np.random.default_rng(seed) for chrom, seed in seeds.items()}
    counts = {rate: 0 for rate in destinations}
    total, variants, header, indices, previous = 0, 0, None, None, {}
    with ExitStack() as stack:
        reader = stack.enter_context(gzip.open(source, 'rt'))
        writers = {rate: stack.enter_context(gzip.open(path, 'wt', compresslevel=1))
            for rate, path in destinations.items()}
        for line in reader:
            if line.startswith('#CHROM'):
                header = line.rstrip().split('\t')
                if len(set(header[9:])) != len(header[9:]):
                    raise ValueError('Duplicate sample IDs')
                indices = sorted(range(len(header) - 9), key=lambda index: header[index + 9])
            if line.startswith('#'):
                for writer in writers.values():
                    writer.write(line)
                continue
            if header is None:
                raise ValueError('Missing genotype header')
            fields = line.rstrip().split('\t')
            chrom, position = fields[0], int(fields[1])
            if chrom not in generators or position <= previous.get(chrom, -1):
                raise ValueError('Unexpected chromosome or unsorted genomic records')
            previous[chrom] = position
            if len(fields) != len(header):
                raise ValueError('Inconsistent sample columns')
            formats = fields[8].split(':')
            genotype_column = formats.index('GT')
            calls = [entry.split(':') for entry in fields[9:]]
            if any(parts[genotype_column] not in {'0|0', '0|1', '1|0', '1|1'} for parts in calls):
                raise ValueError('Expected complete phased biallelic source genotypes')
            uniforms = generators[chrom].random(len(calls))
            for rate, writer in writers.items():
                masked = [parts.copy() for parts in calls]
                for rank, index in enumerate(indices):
                    if uniforms[rank] < rate:
                        masked[index][genotype_column] = '.|.'
                        counts[rate] += 1
                writer.write('\t'.join(fields[:9] + [':'.join(parts) for parts in masked]) + '\n')
            total += len(calls)
            variants += 1
    return {'variants': variants, 'genotypes': total,
        'masked': {str(rate): {'count': count, 'fraction': count / total} for rate, count in counts.items()}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--protocol', type=Path, required=True)
    parser.add_argument('--cohort', type=Path, required=True)
    parser.add_argument('--index', type=int, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):
        parser.error('Full genotype masking requires Slurm')
    protocol = json.loads(args.protocol.read_text())
    item = protocol['cohorts'][args.index]
    args.output.mkdir(parents=True, exist_ok=False)
    source = args.cohort / 'family.normalized.vcf.gz'
    source_hash = sha256(source)
    destinations = {rate: args.output / f'missing_{rate}.vcf.gz' for rate in protocol['missingness_contract']['rates']}
    counts = mask_genotypes(source, destinations,
        {chrom: entry['dropout_contract_only'] for chrom, entry in item['derived_seeds']['chromosomes'].items()})
    checks = []
    for rate, masked in destinations.items():
        prepared = args.output / f'prepared_{rate}.vcf.gz'
        try:
            prepare_vcf(masked, prepared)
        except ValueError as error:
            if 'missing allele' not in str(error) or prepared.exists():
                raise ValueError('Missing-data rejection did not preserve the input contract') from error
            checks.append({'rate': rate, 'rejected_missing_allele': True, 'published_prepared_vcf': False,
                'message': str(error), 'masked_sha256': sha256(masked)})
        else:
            raise ValueError('Missing-data input was unexpectedly accepted')
    if sha256(source) != source_hash:
        raise ValueError('Source genotype input changed')
    save_json(args.output / 'result.json', {'state': 'MISSINGNESS_REJECTION_VERIFIED', 'counts': counts,
        'checks': checks, 'source_sha256': source_hash, 'protocol_sha256': sha256(args.protocol),
        'biological_robustness_qualified': False, 'alleles_imputed': 0, 'sites_deleted': 0})


if __name__ == '__main__':
    main()
