"""Archive explicitly selected frozen assets by content hash, without raw inputs."""
import argparse
import hashlib
import io
import json
from pathlib import Path
import tarfile


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 << 20), b''):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('selection', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    selection = json.loads(args.selection.read_text())
    root = Path(selection['original_root'])
    expected = selection.get('expected_sha256', {})
    records = []
    objects = {}
    for name in selection['files']:
        path = Path(name)
        relative = path.relative_to(root).as_posix()
        digest = sha256(path)
        if name in expected and digest != expected[name]:
            raise ValueError('Frozen binding mismatch: ' + relative)
        records.append({'path': relative, 'sha256': digest,
                        'bytes': path.stat().st_size,
                        'mode': path.stat().st_mode & 0o777})
        objects.setdefault(digest, path)
    manifest = {'schema': 'ghist.frozen-assets.v1',
                'original_root': str(root), 'files': records,
                'purpose': selection['purpose'],
                'known_frozen_bindings_verified': len(expected)}
    manifest_bytes = (json.dumps(manifest, indent=2) + '\n').encode()
    if args.output.exists():
        raise FileExistsError(args.output)
    with tarfile.open(args.output, 'w:gz', compresslevel=4, dereference=True) as archive:
        info = tarfile.TarInfo('manifest.json')
        info.size = len(manifest_bytes)
        info.mode = 0o644
        archive.addfile(info, io.BytesIO(manifest_bytes))
        for digest, path in sorted(objects.items()):
            archive.add(path, arcname='objects/' + digest, recursive=False)
    args.output.with_suffix(args.output.suffix + '.manifest.json').write_bytes(manifest_bytes)
    report = {'archive': args.output.name, 'sha256': sha256(args.output),
              'bytes': args.output.stat().st_size, 'files': len(records),
              'distinct_objects': len(objects),
              'uncompressed_bytes': sum(p.stat().st_size for p in objects.values())}
    args.output.with_suffix(args.output.suffix + '.receipt.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    main()
