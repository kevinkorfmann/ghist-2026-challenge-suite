import hashlib, importlib.metadata as md, json, pathlib, platform, sys
root = pathlib.Path(sys.argv[1])
manifest = pathlib.Path(sys.argv[2]) if len(sys.argv) > 2 else root / 'experiments/readiness_2026_11_20/repair_integration/source_snapshot.json'
obj = json.loads(manifest.read_text())
checks = {p: hashlib.sha256((root/p).read_bytes()).hexdigest() == h for p,h in obj['files_sha256'].items()}
print(json.dumps({'python':sys.version, 'platform':platform.platform(), 'packages':{p:md.version(p) for p in ['numpy','scipy','moments-popgen','dadi','msprime','tskit','scikit-learn','pytest']}, 'source_manifest_sha256':hashlib.sha256(manifest.read_bytes()).hexdigest(), 'source_files_checked':len(checks), 'source_mismatches':[p for p,ok in checks.items() if not ok]},indent=2))
assert all(checks.values())
