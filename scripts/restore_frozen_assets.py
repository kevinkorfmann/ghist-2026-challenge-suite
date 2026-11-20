"""Verify and restore the content-addressed GHIST model/source release bundle."""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import shutil
import tarfile
import tempfile


def restore(archive, destination, verify_only=False):
    destination = Path(destination).resolve()
    with tarfile.open(archive, 'r|gz') as bundle:
        entry = bundle.next()
        if entry is None or entry.name != 'manifest.json' or not entry.isfile() or entry.size > 10_000_000:
            raise ValueError('Invalid manifest')
        manifest = json.load(bundle.extractfile(entry))
        if manifest.get('schema') != 'ghist.frozen-assets.v1':
            raise ValueError('Unsupported archive schema')
        records = manifest['files']
        by_hash = {}
        for record in records:
            path = PurePosixPath(record['path'])
            digest = record['sha256']
            if path.is_absolute() or '..' in path.parts or not path.parts:
                raise ValueError('Unsafe output path')
            if len(digest) != 64 or any(c not in '0123456789abcdef' for c in digest):
                raise ValueError('Invalid object digest')
            target = destination / path
            if not target.resolve().is_relative_to(destination):
                raise ValueError('Output path leaves destination')
            by_hash.setdefault(digest, []).append(record)
        with tempfile.TemporaryDirectory(prefix='ghist-assets-') as temporary:
            scratch = Path(temporary) / 'object'
            seen = set()
            for member in bundle:
                if member.name == 'manifest.json':
                    continue
                digest = member.name.removeprefix('objects/')
                if member.name != 'objects/' + digest or digest not in by_hash or digest in seen:
                    raise ValueError('Unexpected or duplicate object')
                group = by_hash[digest]
                if not member.isfile():
                    raise ValueError('Object is not a regular file')
                seen.add(digest)
                h = hashlib.sha256()
                with bundle.extractfile(member) as source, scratch.open('wb') as output:
                    for block in iter(lambda: source.read(8 << 20), b''):
                        h.update(block)
                        output.write(block)
                if h.hexdigest() != digest or any(scratch.stat().st_size != r['bytes'] for r in group):
                    raise ValueError('Object hash or size mismatch')
                if verify_only:
                    continue
                for record in group:
                    target = destination / record['path']
                    target.parent.mkdir(parents=True, exist_ok=True)
                    if not target.resolve().is_relative_to(destination):
                        raise ValueError('Output path leaves destination')
                    if target.exists():
                        if hashlib.sha256(target.read_bytes()).hexdigest() != digest:
                            raise FileExistsError(f'Refusing to overwrite different file: {target}')
                        continue
                    shutil.copyfile(scratch, target)
                    target.chmod(record.get('mode', 0o644) & 0o777)
            if seen != set(by_hash):
                raise ValueError('Missing archive objects')
        print(f"PASS: {len(records)} paths, {len(by_hash)} unique objects")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', type=Path)
    parser.add_argument('destination', type=Path)
    parser.add_argument('--verify-only', action='store_true')
    args = parser.parse_args()
    restore(args.archive, args.destination, args.verify_only)


if __name__ == '__main__':
    main()
