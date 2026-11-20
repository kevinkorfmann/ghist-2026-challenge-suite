"""Verify corrected VCF export against an existing frozen simulated spectrum."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import tskit

from ghist_inference.checkpoint import atomic_json, load_spectrum_checkpoint
from scripts.validate_bottleneck_readiness import vcf_contract


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--trees', type=Path, required=True)
    parser.add_argument('--spectrum', type=Path, required=True)
    parser.add_argument('--protocol', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    protocol_hash = hashlib.sha256(args.protocol.read_bytes()).hexdigest()
    config = json.loads(args.protocol.read_text())
    blocks, metadata = load_spectrum_checkpoint(args.spectrum, expected={'protocol_sha256': protocol_hash})
    identity = {key: metadata[key] for key in ('case', 'point', 'replicate', 'scenario',
                                              'ancestry_seed', 'mutation_seed', 'optimizer_seed', 'protocol_sha256')}
    ts = tskit.load(args.trees)
    length = config['sequence_length']
    if ts.sequence_length != length or ts.num_samples != 114 or ts.num_individuals != 57:
        raise ValueError('Tree sequence does not match the frozen sample/exposure contract')
    actual_blocks = ts.allele_frequency_spectrum(mode='site', polarised=True, span_normalise=False,
                                                windows=np.linspace(0, length, 51))
    np.testing.assert_array_equal(actual_blocks, blocks)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    result = vcf_contract(ts, length, blocks, args.out.parent, identity, metadata['spectrum_sha256'])
    result['validation_revision'] = 'exact_header_and_bound_provenance_v2'
    result['validation_source_sha256'] = {name: hashlib.sha256(Path(name).read_bytes()).hexdigest() for name in
                                         ('scripts/revalidate_bottleneck_vcf.py', 'scripts/validate_bottleneck_readiness.py',
                                          'src/ghist_bottleneck/sfs_baseline.py', 'src/ghist_inference/checkpoint.py')}
    with args.trees.open('rb') as handle:
        result['tree_sha256'] = hashlib.file_digest(handle, 'sha256').hexdigest()
    atomic_json(args.out, result)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
