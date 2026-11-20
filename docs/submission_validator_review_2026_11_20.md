# Submission validator robustness review — 2026-11-20

The current-source validator now rejects malformed relatedness CSV syntax and
malformed expected-pair domains before reporting verified coverage. The changes
cover local file format and supplied pair coverage only. They do not establish
inference accuracy or prove that a caller supplied the actual input cohort.

## Corrected behavior

- The `--pairs` CSV reader previously accepted duplicate headers and silently
  discarded extra cells. For example, `indID1,indID2,indID1` followed by `x,b,a`
  could certify coverage of `(a,b)` after overwriting the first ID. Headers must
  now be unique and contain both required IDs; row widths must match the header.
  Extra named metadata columns and repeated/reversed expected pairs remain valid.
- Strict CSV parsing rejects unterminated quoted fields and text after a closing
  quote. Previously an unterminated quoted meiosis value could be parsed as an
  integer and accepted. CLI failures produce the usual JSON failure report and
  nonzero exit status.
- Expected pairs must contain two distinct, nonempty string IDs. A string such
  as `"ab"` is no longer interpreted as pair `(a,b)`. Coverage helpers normalize
  surrounding whitespace consistently with the existing submission reader and
  reject IDs that collide after normalization.
- Family coverage validates the actual sample cohort and consumes generators
  once. A valid generator previously lost all its members on the second pass;
  malformed cohorts with duplicate IDs could previously pass a set comparison.
  VCF sample-header validation now uses linear space instead of constructing and
  discarding every pair.
- Supplying expected-pair coverage to Challenges 1–8 now raises a clear error
  instead of silently ignoring it.
- Large finite YAML integers no longer cause `math.isfinite` to raise
  `OverflowError`. Integer finiteness does not require conversion to float.
  This remains a format check, without a new biological upper bound.

Valid relation labels, meiosis mappings, CSV/TSV output schemas, sweep schemas,
and demographic parameter constraints are unchanged. No scientific model,
predictor, challenge configuration, or frozen artifact was edited.

## Verification

The final focused run passed **104 tests**:

```sh
.venv/bin/python -m pytest tests/test_submission_validation.py tests/test_relatedness_submission_contract.py tests/test_task9_release.py tests/test_reliability.py -q
```

Ruff passed on the validator, its CLI, and both assigned contract test files.
Regressions exercise malformed quoting, malformed expected-pair CSVs and JSON
reports, invalid ID domains, generator family partitions, VCF header order,
coverage misuse, large finite YAML integers, and preserved valid CSV forms.

The coordinating agent's combined compatibility run passed **482 tests** across
the submission validator, relatedness contract, downloader, inventory verifier,
Slurm checks, Task 9 release/package, and Task 10/11 complete-pipeline unit tests.
The current Task 9 package tests pass with the explicit package-checksum
argument. This review did not edit the package owner's files.

## Pinned Task 9 snapshot

Before editing, the package manifest at
`experiments/submission_package_2026_11_20/task9/package/package.json` was checked
for validator hashes. The snapshot and manifest were preserved. A read-only
checksum check after the edits confirmed that all **16** snapshot files still
match that manifest. The current-source files intentionally differ:

| File | Pinned snapshot SHA-256 | Current-source SHA-256 |
| --- | --- | --- |
| `src/ghist_submission/validate.py` | `b922bce5f23a9d5ae5d56237fca7e110b664fbc82b735bd4b890c2fc8e8425cc` | `f405f5beaf86c64e1c726b33a97d8cf69aacd10e98c1e9e28b06e2c48f4dee17` |
| `scripts/validate_submission.py` | `747a2041e8e7e00d69678a22ade7b192ee78fa7fe847bf1b3da329d7f35971aa` | `c88e3dbda34254996a070c13432237dc6bebbed7815749dcec1916668128602e` |

`configs/task9_current_release.json` pins relatedness feature-contract source
files, not these validators. Existing operational/scientific evidence remains
tied to its original source snapshot; these changes do not silently re-sign it.
The preserved package still contains its earlier validator and therefore does
not acquire the new malformed-input checks. A package owner can incorporate
the fixes in a separately verified package when appropriate.

No public testing answers, benchmark holdout content, final data, final-data
inference, submissions, publishing, commits, or pushes were used in this review.
