"""Recover documented POS-only correction; never claim current-release byte identity."""
from pathlib import Path
import datetime
import gzip
import hashlib
import json
import re

root = Path('/private/tmp/ghist-final-submissions-20261120/task5')
source = root / 'onesweep.final.original_v1.vcf.gz'
target = root / 'onesweep.final.reconstructed_pos1.vcf.gz'
assert not target.exists()
source_bytes = source.read_bytes()
assert hashlib.md5(source_bytes).hexdigest() == '242445cad013f18cf7c661a869ccd2d1'
assert hashlib.sha256(source_bytes).hexdigest() == '4ffd4b17fde78c182fd79b73ea7ebb246acb077d0bdc75fc04ff65c905377a14'
lengths = {}
count = 0
last = -1
min_pos = None
with gzip.open(source, 'rb') as inp, target.open('xb') as raw_out:
    with gzip.GzipFile(filename='', mode='wb', fileobj=raw_out, mtime=0) as out:
        for line in inp:
            if line.startswith(b'#'):
                if line.startswith(b'##contig='):
                    match = re.fullmatch(rb'##contig=<ID=([^,>]+),length=(\d+)>\r?\n', line)
                    if match:
                        lengths[match[1]] = int(match[2])
                out.write(line)
                continue
            chrom, raw_pos, rest = line.split(b'\t', 2)
            pos = int(raw_pos)
            assert chrom == b'9' and chrom in lengths
            assert 0 <= pos < lengths[chrom] and pos >= last
            last = pos
            if min_pos is None:
                min_pos = pos
            out.write(chrom + b'\t' + str(pos + 1).encode('ascii') + b'\t' + rest)
            count += 1
assert count > 0
# Separate paired pass verifies every original byte except the integer POS field.
verified = 0
with gzip.open(source, 'rb') as original, gzip.open(target, 'rb') as corrected:
    for old in original:
        new = corrected.readline()
        assert new
        if old.startswith(b'#'):
            assert old == new
        else:
            old_fields = old.split(b'\t')
            new_fields = new.split(b'\t')
            assert len(old_fields) == len(new_fields)
            assert old_fields[0] == new_fields[0]
            assert old_fields[2:] == new_fields[2:]
            assert int(new_fields[1]) == int(old_fields[1]) + 1
            verified += 1
    assert not corrected.read(1)
assert verified == count
target_bytes = target.read_bytes()
record = {
    'task': 5,
    'created_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'input_kind': 'reconstructed documented coordinate-only correction of authenticated original final release',
    'source_entity': 'syn75346852', 'source_version': 1,
    'source_path': str(source),
    'source_md5': hashlib.md5(source_bytes).hexdigest(),
    'source_sha256': hashlib.sha256(source_bytes).hexdigest(),
    'target_path': str(target),
    'target_sha256': hashlib.sha256(target_bytes).hexdigest(),
    'target_md5': hashlib.md5(target_bytes).hexdigest(),
    'target_bytes': len(target_bytes),
    'current_release': {'version': 2, 'metadata_md5': 'f8a29f8bb6b05273e81b290f6d0763cd',
                        'release_note': 'Correct VCF POS to 1-based (+1); all other columns unchanged'},
    'download_failure': 'Repeated current-version direct/official ZIP S3 AccessDenied',
    'current_release_bytes_verified': False,
    'current_release_content_equivalence': 'Inferred only from official POS-only release note; exact current bytes unavailable',
    'transform': 'Add exactly 1 to every POS integer; preserve all other VCF column and header bytes; gzip recompressed deterministically',
    'records': count, 'paired_records_verified': verified,
    'header_lengths': {k.decode(): v for k, v in lengths.items()},
    'original_position_range': [min_pos, last],
    'corrected_position_range': [min_pos + 1, last + 1],
    'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'scientific_parameters_changed': False,
    'independent_review_pending': True,
    'external_submission_performed': False,
}
output = root / 'reconstructed_pos1_provenance.json'
assert not output.exists()
output.write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps(record, indent=2))
