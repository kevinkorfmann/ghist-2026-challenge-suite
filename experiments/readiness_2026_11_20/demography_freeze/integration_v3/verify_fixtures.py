import hashlib, json
from pathlib import Path
root = Path.cwd()
manifest = root / 'fixture_snapshot.json'
expected = json.loads(manifest.read_text())
mismatches = [p for p, h in expected.items() if hashlib.sha256((root/p).read_bytes()).hexdigest() != h]
print(json.dumps({'fixture_manifest_sha256': hashlib.sha256(manifest.read_bytes()).hexdigest(), 'fixture_files_checked': len(expected), 'fixture_mismatches': mismatches}, indent=2))
assert not mismatches
