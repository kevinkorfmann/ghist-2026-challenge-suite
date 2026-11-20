"""Snapshot compact synthetic-development evidence without BAM/FASTQ copying.

Run on Betty with --root pointing only at this task's isolated experiment root.
No arbitrary challenge input directories are traversed.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path
import tarfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.root.name != 'lowcoverage-20261120':
        raise ValueError('expected this isolated experiment root')
    prefixes = ('smoke_', 'continuation_', 'iid_gl2_', 'truth_control_', 'gl1_refit_', 'continued_gl1_',
                'full_', 'profile_', 'window_', 'matched_', 'source_v', 'official_templates')
    files = []
    for path in args.root.rglob('*'):
        relative = path.relative_to(args.root)
        if not path.is_file() or '__pycache__' in relative.parts:
            continue
        if len(relative.parts) > 1 and not relative.parts[0].startswith(prefixes):
            continue
        if path.suffix not in ('.py', '.json', '.log', '.sfs', '.arg', '.sbatch', '.yaml', '.txt', '.npy'):
            continue
        if path.stat().st_size > 10_000_000:
            continue
        files.append(path)
    # Retain the explicit release download provenance, without vendoring binaries.
    provenance = args.root/'external/winsfs-release-v0.7.0/provenance.json'
    if provenance.is_file():
        files.append(provenance)
    manifest = {'created_at_utc': datetime.now(timezone.utc).isoformat(), 'root': str(args.root),
                'status': 'development_only_no_submission_readiness', 'files': []}
    with tarfile.open(args.output, 'w:gz') as archive:
        for path in sorted(set(files)):
            content = path.read_bytes()
            relative = str(path.relative_to(args.root))
            manifest['files'].append({'path': relative, 'bytes': len(content), 'sha256': hashlib.sha256(content).hexdigest()})
            info = tarfile.TarInfo(relative)
            info.size = len(content)
            archive.addfile(info, io.BytesIO(content))
        content = (json.dumps(manifest, indent=2)+'\n').encode()
        info = tarfile.TarInfo('artifact_manifest.json')
        info.size = len(content)
        archive.addfile(info, io.BytesIO(content))
    print(json.dumps({'archive': str(args.output), 'sha256': hashlib.sha256(args.output.read_bytes()).hexdigest(),
                      'files': len(manifest['files'])}))


if __name__ == '__main__':
    main()
