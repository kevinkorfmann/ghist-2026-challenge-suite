# Download and inventory reliability review — November 20, 2026

This review covers `scripts/download_synapse.py` and
`scripts/verify_data_inventory.py`. Verification used synthetic temporary files
and a mocked Synapse client. No authenticated download, existing data-file audit,
final-data execution, challenge-answer read, manifest replacement, submission,
commit, or push was performed.

## Repairs and operating behavior

Downloads previously wrote directly into the version directory and replaced its
manifest on every invocation. They now use a unique staging directory, validate
the returned entity ID, version, and file handle against the resolved metadata,
and require a regular downloaded file inside that staging directory. When
Synapse provides size or MD5 metadata, the downloaded bytes must match it. An
optional `--expected-sha256` accepts a separately trusted SHA-256 digest; a
mismatch prevents publication. The recorded `sha256` by itself remains an
observation of local bytes, not independent proof of their expected contents.

An existing version directory is rejected without altering it, including when a
second invocation creates that directory while the first is downloading. The
new directory is created exclusively after verification, and its manifest is
published last. A recoverable publication failure removes only the new directory
created by that invocation. A process kill or power loss can still leave a
staging directory or an incomplete version directory without a manifest. Such a
version directory is rejected on retry; inspect it and use a different
`--data-root` if another download is needed. Automatic overwrite or resume is
intentionally absent. Publication does not provide a power-loss durability
guarantee.

Registry paths must remain relative beneath `challenges/`. Downloaded files
named `manifest.json`, missing download paths, paths outside staging, and symlink
payloads are rejected. External `file://` entities that return an outside path
therefore require a separate intake workflow. Available Synapse metadata is
checked; an absent size or MD5 does not manufacture an expected checksum.

Inventory verification previously considered a missing file with a null expected
digest to be a match (`None == None`), treated an empty inventory as success, and
stopped at the first read error. It now requires a valid SHA-256 string and a
nonempty inventory, reports individual malformed/missing/unreadable/mismatched
entries, and continues to later files. Hexadecimal digests are case insensitive.
`exists: null` means a malformed entry was rejected before its file was checked.
Malformed YAML fields and NUL-containing paths remain reportable as JSON. A
failed audit exits with status 1 after emitting its report.

Report files are replaced atomically from a temporary file, and a report path
that resolves to the inventory or any inventoried input is rejected. An
unwritable or invalid report destination still raises an error and exits
nonzero; it is not reported as a completed audit. Inventory input paths retain
their existing interpretation relative to the invocation's working directory.
Matching the historical inventory still does not establish that a file is the
latest Synapse revision.

## Verification

The following commands passed using Python 3.12.11, pytest 9.2.0, Ruff 0.17.6,
PyYAML 6.0.3, and synapseclient 4.15.0:

```bash
.venv/bin/python -m pytest -q tests/test_download_synapse.py tests/test_data_inventory.py
# 39 passed
.venv/bin/ruff check scripts/download_synapse.py scripts/verify_data_inventory.py \
  tests/test_download_synapse.py tests/test_data_inventory.py
# All checks passed!
git diff --check -- scripts/download_synapse.py scripts/verify_data_inventory.py \
  tests/test_download_synapse.py
```

Regression coverage includes interrupted transfers, existing and concurrently
created destinations, publication failure, identity/version/file-handle mismatch,
size/MD5/SHA-256 mismatch, escaped paths, missing files, malformed inventories,
per-file I/O errors, complete failure reports, and attempted report overwrite of
the input or inventory. Installed Synapse SDK source was inspected, and a real
`synapseclient.File` was constructed offline: `id` and `dataFileHandleId` are
strings, `versionNumber` and the `fileSize` alias are integers, and the `md5`
alias is a string. These agree with the new checks. Live service behavior,
authentication, and scientific inference were not validated by these tests.

## Validation provenance

Before changes, the three existing source/test files exactly matched their
entries in `experiments/reliability_2026_11_20/verification.json`. That historical
record is preserved and no longer describes these three current files. Its
other results have not been regenerated or represented as validation of this
repair. No scientific inference implementation or freeze record was changed.
The new inventory test was not present in that prior source inventory.

Current SHA-256 values for the files exercised above are:

```json
{
  "scripts/download_synapse.py": "5c5776de2818a5d718b4e1552fcf4ccad7f442928d6614ca09993a9e33fcdfe4",
  "scripts/verify_data_inventory.py": "b22c54f4f831a309ec7cc65547a70059c36fd38305c1546605810a585335ad19",
  "tests/test_download_synapse.py": "826e974974e7aac015f490b406707d33afe7a84a50aebafec2a3fd0051a7878c",
  "tests/test_data_inventory.py": "d1c8edb594fff53a8fe618b56418f36d6f363b7540135b1b92265a2237170527"
}
```
