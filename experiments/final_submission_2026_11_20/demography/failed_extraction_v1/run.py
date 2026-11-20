"""Final-data best-effort inference; never freezes methods or submits externally."""
from __future__ import annotations
import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import sys
import time
import traceback

ROOT = Path(__file__).resolve().parent
SNAPSHOT = ROOT / 'source_snapshot'
sys.path[:0] = [str(SNAPSHOT / 'src'), str(SNAPSHOT)]
import numpy as np
import pysam
from ghist_inference.checkpoint import atomic_json, spectrum_digest


def sha(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            result.update(chunk)
    return result.hexdigest()


def verify_source(expected_lock_sha256=None):
    lock_hash = sha(ROOT / 'method_lock.json')
    if expected_lock_sha256 is not None and lock_hash != expected_lock_sha256:
        raise ValueError('Method lock changed during execution')
    manifest = json.loads((ROOT / 'source_manifest.json').read_text())
    for name, digest in manifest.items():
        if sha(SNAPSHOT / name) != digest:
            raise ValueError(f'Source changed: {name}')
    protocol = json.loads((ROOT / 'method_lock.json').read_text())
    if sha(ROOT / 'method_lock.json') != lock_hash:
        raise ValueError('Method lock changed during parsing')
    if sha(ROOT / 'source_manifest.json') != protocol['source_manifest_sha256']:
        raise ValueError('Source manifest differs from lock')
    if sha(Path(__file__)) != protocol['runner_sha256']:
        raise ValueError('Runner differs from lock')
    return protocol


def input_record(task):
    record = json.loads((ROOT / 'inputs' / f'task{task}.json').read_text())
    for item in record['files'].values():
        path = Path(item['path'])
        if not path.is_absolute():
            path = ROOT / path
        if sha(path) != item['sha256'] or path.stat().st_size != item['bytes']:
            raise ValueError('Final input identity differs from authenticated handoff')
        item['resolved_path'] = str(path)
    return record


def vcf_preflight(path, task):
    with pysam.VariantFile(str(path)) as source:
        lengths = {name: source.header.contigs[name].length for name in source.header.contigs}
        samples = list(source.header.samples)
        if not samples or len(set(samples)) != len(samples):
            raise ValueError('Samples absent or duplicated')
        if not lengths or any(v is None or v <= 0 for v in lengths.values()):
            raise ValueError('Contig lengths absent or invalid')
        ranks = {name: i for i, name in enumerate(lengths)}
        previous = None
        records = biallelic = invalid_ancestral = 0
        for row in source:
            point = (ranks[row.contig], row.pos)
            if not 1 <= row.pos <= lengths[row.contig] or (previous is not None and point <= previous):
                raise ValueError('Out-of-range, unsorted or duplicate VCF locus')
            previous = point
            records += 1
            if len(row.alleles or ()) != 2 or any(a not in {'A','C','G','T'} for a in row.alleles):
                continue
            biallelic += 1
            if task != 1:
                aa = row.info.get('AA')
                if isinstance(aa, tuple): aa = aa[0] if len(aa) == 1 else None
                invalid_ancestral += aa not in row.alleles
        if not biallelic or invalid_ancestral:
            raise ValueError('Empty SNP set or invalid ancestral states')
    return {'contig_lengths': lengths, 'diploid_samples': len(samples), 'records': records,
            'biallelic_nucleotide_records': biallelic, 'invalid_ancestral_records': invalid_ancestral,
            'sorted_unique_loci': True, 'sequence_length': sum(lengths.values())}


def extract(task):
    lock_hash = sha(ROOT / 'method_lock.json')
    verify_source(lock_hash)
    inputs = input_record(task)
    out = ROOT / 'results' / f'task{task}'
    out.mkdir(parents=True, exist_ok=True)
    if (out / 'spectrum.json').exists() or (out / 'spectrum.npz').exists():
        raise ValueError('Refusing to overwrite existing extraction')
    vcf = Path(inputs['files']['vcf']['resolved_path'])
    audit = vcf_preflight(vcf, task)
    length = audit['sequence_length']
    if task == 1:
        from ghist_bottleneck.sfs_baseline import projected_sfs_from_vcf, spectrum_exposure
        spectrum, extraction = projected_sfs_from_vcf(vcf, polarized=False)
        length = spectrum_exposure(extraction, length, None)
    elif task == 2:
        from ghist_secondary_contact.vcf_fast import scan_secondary_contact_vcf
        if len(audit['contig_lengths']) != 1: raise ValueError('Task2 requires one complete contig')
        spectrum, extraction = scan_secondary_contact_vcf(vcf, sequence_length=length)
    elif task == 3:
        from ghist_growth_bgs.audit import audit_growth_resources
        from ghist_growth_bgs.spectrum import spectrum_from_counts
        extraction = audit_growth_resources(vcf, Path(inputs['files']['coding']['resolved_path']), Path(inputs['files']['map']['resolved_path']))
        labels = ('noncoding_1kb_10kb', 'noncoding_10kb_100kb')
        spectrum = sum(spectrum_from_counts(extraction['distance_strata'][x]['unfolded_sfs'], extraction['expected_haplotypes']) for x in labels)
        length = sum(extraction['physical_bp_by_distance_stratum'][x] for x in labels)
    else:
        from ghist_split_migration_bgs.spectrum import stratified_joint_spectra_from_vcf
        from ghist_growth_bgs.audit import load_coding_intervals, load_recombination_map, physical_stratum_lengths
        coding = Path(inputs['files']['coding']['resolved_path'])
        recombination = Path(inputs['files']['map']['resolved_path'])
        spectra, extraction = stratified_joint_spectra_from_vcf(vcf, coding)
        lengths, _, _ = physical_stratum_lengths(audit['contig_lengths'], load_coding_intervals(coding), load_recombination_map(recombination), lambda _: 'unused')
        labels = ('noncoding_1kb_10kb', 'noncoding_10kb_100kb')
        spectrum = sum(spectra[x] for x in labels)
        length = sum(lengths[x] for x in labels)
        extraction['physical_bp_by_distance_stratum'] = lengths
    after = input_record(task)
    if after != inputs: raise ValueError('Input record changed during extraction')
    if length <= 0 or not np.isfinite(spectrum).all() or spectrum.sum() <= 0:
        raise ValueError('Invalid spectrum or exposure')
    verify_source(lock_hash)
    target = out / 'spectrum.npz'
    with (out / 'spectrum.tmp').open('wb') as stream:
        np.savez_compressed(stream, spectrum=spectrum, sequence_length=float(length))
    os.replace(out / 'spectrum.tmp', target)
    atomic_json(out / 'spectrum.json', {'task': task, 'input_record': inputs, 'audit': audit,
                'extraction': extraction, 'spectrum_sha256': spectrum_digest(spectrum), 'archive_sha256': sha(target),
                'sequence_length': float(length), 'method_lock_sha256': lock_hash,
                'source_manifest_sha256': sha(ROOT / 'source_manifest.json')})


def fit(index):
    lock_hash = sha(ROOT / 'method_lock.json')
    protocol = verify_source(lock_hash)
    case = protocol['fits'][index]
    task = case['task']
    inputs = input_record(task)
    out = ROOT / 'results' / f'task{task}'
    record_path = out / (case['id'] + '.json')
    if record_path.exists(): raise ValueError('Refusing to overwrite a fit')
    spectrum_record = json.loads((out / 'spectrum.json').read_text())
    if spectrum_record['input_record'] != inputs or spectrum_record['method_lock_sha256'] != lock_hash:
        raise ValueError('Spectrum is not linked to current locked inputs/method')
    archive_path = out / 'spectrum.npz'
    if sha(archive_path) != spectrum_record['archive_sha256']: raise ValueError('Spectrum archive changed')
    with np.load(archive_path) as archive:
        spectrum, length = archive['spectrum'].copy(), float(archive['sequence_length'])
    if spectrum_digest(spectrum) != spectrum_record['spectrum_sha256']: raise ValueError('Spectrum digest mismatch')
    start = time.monotonic()
    result = {'state': 'RUNNING', 'case': case, 'slurm_job_id': os.environ.get('SLURM_JOB_ID'),
              'method_lock_sha256': lock_hash, 'spectrum_sha256': spectrum_record['spectrum_sha256'],
              'input_vcf_sha256': inputs['files']['vcf']['sha256'], 'scientifically_validated': False, 'submission_count': 0,
              'software': {p: importlib.metadata.version(p) for p in ('numpy','scipy','moments-popgen','dadi','pysam')}}
    atomic_json(record_path, result)
    try:
        args = dict(case['settings'])
        if 'projection' in args and isinstance(args['projection'], list): args['projection'] = tuple(args['projection'])
        if task == 1:
            if case['engine'] == 'moments':
                from ghist_bottleneck.sfs_baseline import fit_two_epoch_spectrum as function
            else:
                from ghist_bottleneck.dadi_reference import fit_two_epoch_dadi_spectrum as function
        elif task == 2:
            from ghist_secondary_contact.spectrum import fit_secondary_contact_spectrum as function
        elif task == 3:
            from ghist_growth_bgs.spectrum import fit_growth_spectrum as function
        else:
            from ghist_split_migration_bgs.spectrum import fit_isolation_with_migration_spectrum as function
        value = function(spectrum, sequence_length=length, **args)
        value['sequence_length'] = length
        value['mutation_rate'] = args['mutation_rate']
        result.update(state='COMPLETED', fit=value)
        atomic_json(record_path, result)  # Preserve finite result before subsequent attestation.
        if input_record(task) != inputs: raise ValueError('Input changed during fitting')
        verify_source(lock_hash)
        if sha(archive_path) != spectrum_record['archive_sha256']: raise ValueError('Spectrum changed during fitting')
        result['postfit_identity_verified'] = True
    except Exception as error:
        result.update(state='FAILED', failure=f'{type(error).__name__}: {error}', traceback=traceback.format_exc())
        raise
    finally:
        result['elapsed_seconds'] = time.monotonic() - start
        atomic_json(record_path, result)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=('extract','fit'))
    parser.add_argument('index', type=int)
    args = parser.parse_args()
    if args.action == 'extract': extract(args.index)
    else: fit(args.index)
