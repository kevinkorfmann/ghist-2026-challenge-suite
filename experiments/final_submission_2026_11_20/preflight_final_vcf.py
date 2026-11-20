"""Read-only VCF structure audit; no prediction, fitting, or truth access."""
import argparse
import collections
import gzip
import hashlib
import json
from datetime import UTC, datetime
from pathlib import Path


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def audit(source, output, expected):
    assert not output.exists()
    before = digest(source)
    assert before == expected
    samples = None
    chromosomes = {}
    genotype_counts = collections.Counter()
    closed = set()
    previous = None
    total = 0
    with gzip.open(source, 'rt') as handle:
        for line in handle:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM\t'):
                assert samples is None
                samples = line.rstrip('\n').split('\t')[9:]
                assert samples and len(samples) == len(set(samples))
                continue
            assert samples is not None and not line.startswith('#')
            fields = line.rstrip('\n').split('\t')
            assert len(fields) == 9 + len(samples)
            chrom, position = fields[0], int(fields[1])
            assert position > 0
            if previous and chrom != previous:
                closed.add(previous)
            assert chrom not in closed
            previous = chrom
            row = chromosomes.setdefault(chrom, {
                'sites': 0, 'complete_sites': 0, 'missing_genotypes': 0,
                'unphased_called_genotypes': 0, 'multiallelic_sites': 0,
                'min_position': position, 'max_position': 0,
            })
            assert position > row['max_position'], (chrom, position)
            row['max_position'] = position
            row['sites'] += 1
            row['multiallelic_sites'] += ',' in fields[4]
            formats = fields[8].split(':')
            assert 'GT' in formats
            gt_index = formats.index('GT')
            missing = 0
            for sample in fields[9:]:
                values = sample.split(':')
                assert gt_index < len(values)
                gt = values[gt_index]
                genotype_counts[gt] += 1
                alleles = gt.replace('|', '/').split('/')
                assert len(alleles) == 2
                assert all(a == '.' or a.isdecimal() for a in alleles)
                if '.' in alleles:
                    missing += 1
                elif '|' not in gt:
                    row['unphased_called_genotypes'] += 1
            row['missing_genotypes'] += missing
            row['complete_sites'] += missing == 0
            total += 1
    assert total and samples
    assert digest(source) == before
    result = {
        'verified_at_utc': datetime.now(UTC).isoformat(),
        'scope': 'Final-input structure only; no model, final truth, or label inference.',
        'input_path': str(source), 'input_sha256': before,
        'script_sha256': digest(Path(__file__)),
        'sample_count': len(samples), 'sample_ids': samples, 'sites': total,
        'chromosomes': chromosomes, 'genotype_counts': dict(genotype_counts),
        'source_unchanged_after_read': True,
    }
    output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'sample_count': len(samples), 'sites': total,
                      'chromosomes': len(chromosomes),
                      'missing_genotypes': sum(x['missing_genotypes'] for x in chromosomes.values()),
                      'unphased_called_genotypes': sum(x['unphased_called_genotypes'] for x in chromosomes.values())}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--expected-sha256', required=True)
    args = parser.parse_args()
    audit(args.input, args.output, args.expected_sha256)
