"""Export already completed synthetic operational evidence, without model selection."""
import hashlib
import json
import sys
import tarfile
from pathlib import Path


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


root = Path(sys.argv[1]).resolve()
smoke = root / 'pragmatic_corrected_context_smoke_v1'
audit = json.loads((smoke / 'smoke_audit.json').read_text())
if not audit['family_partition_unchanged'] or not all(v['passed'] for v in audit['comparisons'].values()):
    raise ValueError('Synthetic operational parity did not pass')
for relative, digest in audit['output_sha256'].items():
    if sha(smoke / relative) != digest:
        raise ValueError(f'Changed smoke evidence: {relative}')
files = [root / name for name in (
    'pragmatic_corrected_base_candidate_v1.json', 'pragmatic_corrected_context_candidate_v1.json',
    'pragmatic_corrected_bindings_v1.py', 'prepare_pragmatic_corrected_v1.py',
    'smoke_pragmatic_corrected_v1.py', 'smoke_pragmatic_corrected_v1.sbatch',
    'pragmatic_corrected_operational_source_v1.json', 'test_pragmatic_corrected_v1.py',
    'freeze_and_smoke_pragmatic_corrected_v1.sbatch', 'pragmatic_corrected_inference_v1.sbatch',
    'export_pragmatic_corrected_operational_v1.py',
    'pragmatic_corrected_context_smoke_v1/smoke_audit.json',
    'pragmatic_corrected_context_smoke_v1/alias_and_order_map.json')]
files.extend(p for p in sorted((smoke / 'adapter').iterdir()) if p.suffix in ('.json', '.csv'))
manifest = {'schema': 'task11-corrected-operational-export-v1',
            'scope': 'Completed synthetic alias/order parity; no final candidate selection or challenge data',
            'sha256': {str(p.relative_to(root)): sha(p) for p in files}}
manifest_path = root / 'pragmatic_corrected_operational_export_v1.json'
with manifest_path.open('x') as stream:
    stream.write(json.dumps(manifest, indent=2) + '\n')
files.append(manifest_path)
archive = root / 'pragmatic_corrected_operational_bundle_v1.tar.gz'
with tarfile.open(archive, 'x:gz') as tar:
    for path in files:
        tar.add(path, arcname=str(path.relative_to(root)), recursive=False)
print(json.dumps({'archive': str(archive), 'sha256': sha(archive), 'files': len(files)}))
