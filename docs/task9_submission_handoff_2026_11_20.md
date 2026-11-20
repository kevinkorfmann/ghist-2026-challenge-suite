# Task 9 submission handoff — 2026-11-20

**The model is validated; a current final-data prediction is still pending.**
The current handoff is the versioned **v2** offline source/evidence package,
with canonical `entity_id` provenance support and the current shared validator.
It does not download or read current final data, run current final inference,
upload a file, or submit anything.

The package is at
`experiments/submission_package_2026_11_20/task9_v2/package/`. Its `package.json`
records SHA-256 values for 16 files, release
`2026-11-20-ten-labels`, 60,266 training pairs and 20,023 independent-confirmation
pairs. It contains prediction/validation source, release/configuration and
compact evidence; large models and maps remain on Betty. Final input version,
checksum, path and output identity are `null`, not inferred from historical
artifacts. A [methods draft](task9_methods_draft_2026_11_20.md) is ready for those
fields to be completed after the actual run.

For transfer,
`experiments/submission_package_2026_11_20/task9_v2/task9-package-v2.tar.gz`
contains the same snapshot (23,459 bytes; SHA-256
`f1df68f7e5f2397fd5b7d2540d330c601034ea66f54d199e93b7a2c3c95df0de`).
The package manifest SHA-256 is
`4c0ec8f9ae12b00b77506b803812c3378f1d43bf873ff206a1db2bffeca5da9e`.
The adjacent `SHA256SUMS` pins the archive, manifest and verification records.
Extracting the archive creates `task9-package-v2/`. This exact archive is staged
on Betty at
`/vast/projects/smathi/cohort/korfmann/ghist2026/reviews/task9-offline-package-v2-20261120/`.

The original `task9/` package and all its 21 files are preserved byte-for-byte.
Its archive SHA-256 remains
`5f1d1bd71c2327de4af7ef17a50a09b287445656d66d51f1eac14d62d8ade697`
and manifest SHA-256 remains
`017f91c139aeaa3e013cd4147e3e846b0ada14c6407497952114e6c75cc993b1`.
Those historical files predate the current provenance/validator compatibility
changes and are not the package selected by this handoff. Their original test
and read-only verification records retain their original meaning.

Check the archive against the checksum in this trusted handoff **before extracting
or executing its Python code**. A checksum read only from the transferred archive
or its adjacent `SHA256SUMS` is not an independent trust anchor. Package `verify`
without an externally supplied checksum establishes internal consistency only;
the `preflight` and `run` CLIs require the trusted manifest checksum explicitly.
These are reproducibility/integrity checks, not protection against malicious
package code: the CLI begins executing before its own verification runs. They
are not digital signatures and assume a trusted checkout/handoff and
an output directory whose contents are not being maliciously replaced during a run.

The v2 read-only Betty preflight **passed** in Slurm job **10851132**:
18 seconds, 2 CPUs, 1,598,148 KiB peak RSS, exit `0:0`. It verified these pinned
files and runtime, including loading the model to check its class/feature
contract. A separate post-preflight package inventory check also passed.

| Item | Verified identity |
|---|---|
| Model | `artifacts/task9_releases/2026-11-20-ten-labels/model.joblib`, 442,225,843 bytes |
| Model SHA-256 | `7fc54aa5f0dff1d8dc01364cf781d7a750970e1558a525c2fbdebcbc04a0611c` |
| Map | `external/ped-sim-map/refined_mf.simmap`, 37,825,969 bytes |
| Map SHA-256 | `e5cf4af7f914f9d675c3e3c5d7648775aefd456710dcab826195148ee7aadae1` |
| Runtime | CPython 3.11.11, Linux x86_64, all six release dependency versions match |
| Interpreter | `/vast/projects/smathi/cohort/korfmann/ghist2026/.venv/bin/python` |

The exact dependency versions are in
`requirements/task9-release-2026-11-20.txt`. The macOS workspace environment is
Python 3.12.11 and is suitable for offline packaging/tests, but it is not the
validated inference runtime. The future preflight refuses that runtime and
checks all pinned dependencies, release source, model/map checksums, model
classes and feature order before inference. The predictor parses the map later;
preflight itself verifies map content identity rather than parsing map coordinates.

Offline preparation and verification are reproducible from the checkout:

```bash
python3 scripts/prepare_task9_submission_package.py prepare \
  --checkout . --release configs/task9_current_release.json \
  --out /absolute/path/to/new-task9-package
python3 /absolute/path/to/new-task9-package/scripts/prepare_task9_submission_package.py \
  verify --package /absolute/path/to/new-task9-package
```

Both package and run destinations must be new directories with no symlink path
components; use their canonical absolute parent paths. The tool rejects
historical freeze manifests, default/unversioned model paths, source/evidence
checksum changes, inconsistent model/runtime/evidence metadata, source symlinks,
unexpected files including unpinned bytecode, and reused output directories.
Bytecode writes are disabled while the package runs. It never searches for a
candidate CSV or takes a previously generated candidate as input.

For the later final-data action, use the following commands from the recorded
Linux runtime **inside an appropriate Slurm allocation**. The v2 package is
already staged at the path above. `PACKAGE`, `INPUT_CSV`,
`INPUT_PROVENANCE_JSON` and `RUN_DIR` must be explicit absolute paths. The final
CSV and provenance paths are intentionally pending until the deferred final-data
step takes place; no pre-existing final input is assumed here.

```bash
source /etc/profile
module load python/3.11.11
TASK_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026
PYTHON="$TASK_ROOT/.venv/bin/python"
SEX_MAP="$TASK_ROOT/external/ped-sim-map/refined_mf.simmap"
PACKAGE="$TASK_ROOT/reviews/task9-offline-package-v2-20261120/task9-package-v2"
PACKAGE_SHA256=4c0ec8f9ae12b00b77506b803812c3378f1d43bf873ff206a1db2bffeca5da9e
"$PYTHON" "$PACKAGE/scripts/prepare_task9_submission_package.py" preflight \
  --package "$PACKAGE" --artifact-root "$TASK_ROOT" --sex-map "$SEX_MAP" \
  --expected-package-sha256 "$PACKAGE_SHA256"

# This separate command reads the supplied final input and performs inference.
: "${INPUT_CSV:?Set the verified current final segment CSV}"
: "${INPUT_PROVENANCE_JSON:?Set the recorded Synapse entity/version/checksum JSON}"
: "${RUN_DIR:?Set a new absolute output directory}"
"$PYTHON" "$PACKAGE/scripts/prepare_task9_submission_package.py" run \
  --package "$PACKAGE" --artifact-root "$TASK_ROOT" --sex-map "$SEX_MAP" \
  --observed "$INPUT_CSV" --input-provenance "$INPUT_PROVENANCE_JSON" \
  --run-dir "$RUN_DIR" --workers 2 --expected-package-sha256 "$PACKAGE_SHA256"
```

The provenance JSON must contain `entity_id: "syn75202631"` (as emitted by the
downloader; legacy `entity` is also accepted), the actual positive
integer `version`, and the SHA-256 of that downloaded version. The version must
come from Synapse's file record; a checksum alone cannot establish an entity's
version or whether that version is current. The offline wrapper does not query
Synapse: verify the current entity/version separately when obtaining the final
input. No placeholder provenance file is supplied that could be mistaken for
a verified record.

The `run` command preserves a snapshot as `input.csv` and
`input_provenance.json`, derives `expected_pairs.csv` from validated raw segments,
and invokes the frozen predictor using `--release` and `--dataset final`.
`commands.json` records its exact command arguments. The standalone validator
then checks the resulting file against those expected pairs. Only after identity
and coverage checks succeed, with the input snapshot, roster, model/map hashes and
package identity still matching their pre-run records, is `handoff.json` written with status
`LOCAL_FINAL_FILE_VALIDATED_NOT_SUBMITTED`, output checksum and checksums for all
existing run files except the handoff report itself. A failure leaves diagnostics in its fresh directory; use a new
directory after fixing the cause.

The future prediction file is exactly
`$RUN_DIR/prediction/task9_final_candidate.csv`. Keep `handoff.json`,
`validation.json`, `preflight.json`, predictor/validator logs and the input
provenance as its evidence. The local validator checks format and observed-pair
coverage; it does not establish final-data accuracy or authenticate Synapse
metadata. Before actual submission, complete the pending fields in the methods
draft and use this CSV/checksum rather than a historical candidate. Upload and
submission remain separate actions and are not implemented by this tool.

Validation of the v2 handoff: **67 tests passed** (50 package tests plus 17
existing release tests), including a synthetic ten-class model, generated map/segments,
full predictor invocation, separate validator invocation, provenance failures,
source/metadata/evidence tampering, runtime mismatch, dangling symlinks and stale
candidate refusal. End-to-end failures cover duplicate/reversed-duplicate
segments, HBD, nonfinite coordinates, inconsistent lengths and missing map
chromosomes. Mutations of model/map accounting, input, roster and package after
prediction are rejected without writing a successful handoff. A characterization
test confirms the frozen feature-guard limits described in the methods draft;
the review does not change those validated features. All tests use temporary
synthetic fixtures for inference. Ruff
checks passed. The current-source compatibility rerun is copied into `task9_v2/`
as `compatibility.json` and `compatibility_operational.json`: job 10850916 reran
four production cases on the same 20,023 synthetic pairs, with every output CSV
byte-identical to its frozen counterpart. This is compatibility evidence, not
new independent scientific confirmation. Every applicable predictor/validator
source hash matches the v2 package.

The separate v2 preflight used the staged package, verified all runtime versions
and model/map hashes, and loaded the pinned model to check its labels and feature
order. It read no challenge input, generated no predictions, copied no 442 MB
model locally and submitted nothing. Its Slurm script, output, accounting and
[verification record](../experiments/submission_package_2026_11_20/task9_v2/local_verification.json)
are stored separately under `task9_v2/`; the original package's read-only hash
check is not represented as this new preflight.
