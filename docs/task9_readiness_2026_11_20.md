# Task 9 operational release candidate

Task 9 has a separately versioned ten-label release candidate. It preserves the
historical forest settings, uses only the original calibration simulations, and
repairs the retired-class probability-support failure. It is selected explicitly
with `--release`; historical model artifacts and freezes remain unchanged. No
current challenge-data predictions or submissions are part of this validation.

## Input domain

The [official relatedness specification](https://www.synapse.org/Synapse:syn74803798/wiki/640950)
describes exact Ped-sim IBD1/IBD2 segments, without HBD or ancestral sharing beyond
the pedigree, using sex-specific maps and crossover interference. Task 9 therefore
does not need the genotype generation, phasing and IBD detector validation that
Tasks 10–12 require. The
[submission tutorial](https://www.synapse.org/Synapse:syn74803798/wiki/640752)
requires individual IDs and a relationship and/or meiosis count.

The supplied data are a segment CSV. The specification lists no separate sample
or pair roster for pairs sharing no IBD. We interpret the operational domain as
every distinct unordered pair present in the CSV; this is an inference from the
available input contract, not an explicit statement about zero-IBD pairs in the
wiki. The predictor validates coverage against the raw input before feature
extraction. An optional independent pair list must agree exactly; extra unlisted
pairs cause an explicit error. It does not label all combinations of IDs or
invent identities absent from the input. The scientific accuracy estimates
remain conditional on pairs with observed IBD.

## Evidence and selection

The corrected model was fitted on 60,266 original calibration pairs using the
historical random forest hyperparameters, with the ten current labels. Neither
public testing answers nor new evaluation outcomes selected hyperparameters.
The two untouched confirmation datasets were generated after this fitted model
was fixed. Their scores are:

| Independent confirmation | Pairs | Exact accuracy | Balanced accuracy | Meioses accuracy |
|---|---:|---:|---:|---:|
| Matched interference, seed 2026112307 | 10,038 | 89.10% | 83.35% | 89.41% |
| Poisson crossover shift, seed 2026112308 | 9,985 | 87.20% | 81.48% | 87.71% |

The preserved historical gates require at least 10,000 matched evaluation pairs,
80% exact accuracy, 75% balanced accuracy and 84% meioses accuracy. The numerical
accuracy gates also pass on the shifted confirmation. All map-length stresses
retain valid accepted-class probability support. These gates are minimum release
criteria; they do not make every relationship inference correct.

`validate_task9_operationally.py` anonymizes simulation IDs, creates official CSV
inputs and independent expected-pair lists, invokes the real production predictor,
then invokes the submission-validator CLI separately. It also reverses pair order,
reorders input rows and columns, and adds a UTF-8 BOM. Predictions must agree across
these equivalent inputs and reproduce the fixed scientific confirmation scores.
The versioned release path repeats these checks after the model is staged.

Betty job 10850210 completed successfully in 3 minutes 58 seconds with four CPUs.
All eight production executions passed, covering 20,023 distinct confirmation
pairs across two datasets. Both model-selection paths reproduced the original
scores exactly and every pair was present once in each output. The 54 targeted
local tests, critical Ruff checks and Slurm shell syntax checks also passed.

The later compatibility run, Betty job 10850916, repeated four versioned
production executions using the current manifest-compatible predictor and
shared validator. Every prediction CSV is byte-identical to its original
confirmation output. This reuses the same 20,023 pairs and is not new independent
biological evidence. The separately rebuilt
[v2 handoff package](task9_submission_handoff_2026_11_20.md) passed its pinned
runtime/model/map preflight in job 10851132; the original package is preserved.

The manifest pins the model checksum, feature order, feature-extractor and label
code, genetic map, genome normalization and scikit-learn version. A mismatch
fails before prediction. The predictor rejects duplicate segments, malformed or
nonfinite coordinates, unsupported HBD, invalid maps and missing pairs. Validated
output replaces a prior output only after the new file passes validation. Final
challenge input requires a matching Synapse entity, version and checksum record.

## Reproduction

The job script and compact results are under
[`experiments/readiness_2026_11_20/task9`](../experiments/readiness_2026_11_20/task9/).
Large models and synthetic CSVs remain on Betty. Run the recorded Slurm script
from the isolated readiness snapshot with the original scikit-learn 1.9.0
environment. The manifest is `configs/task9_current_release.json`; model paths
resolve relative to the explicitly supplied artifact root.

The general repository lock uses a different scikit-learn version. Restore the
validated runtime separately using the exact recorded
[release requirements](../requirements/task9-release-2026-11-20.txt). From a clean
checkout on Linux x86_64 with `uv` installed:

```bash
uv python install 3.11.11
uv venv --python 3.11.11 .venv-task9-release
uv pip install --python .venv-task9-release/bin/python \
  -r requirements/task9-release-2026-11-20.txt
.venv-task9-release/bin/python -c \
  'import platform, sklearn; assert platform.python_version() == "3.11.11"; assert sklearn.__version__ == "1.9.0"'
```

Invoke this interpreter directly rather than `uv run`, which can synchronize the
general project lock. On Betty, `TASK_ROOT` is
`/vast/projects/smathi/cohort/korfmann/ghist2026`; the recorded model remains at
`$TASK_ROOT/artifacts/task9_releases/2026-11-20-ten-labels/model.joblib`. For another
machine, copy that model and the recorded genetic map to the corresponding
artifact-root location; the release loader verifies their SHA-256 values. This
requirements file reconstructs the recorded dependency versions; a clean install
on another platform has not been scored independently.

For a future authorized current-data run, select the release explicitly:

```bash
PYTHONPATH=src:. .venv-task9-release/bin/python scripts/predict_task9_current.py \
  --release configs/task9_current_release.json --artifact-root "$TASK_ROOT" \
  --observed "$INPUT_CSV" --sex-map "$SEX_MAP" --dataset final \
  --input-provenance "$INPUT_PROVENANCE_JSON" --workers 2 --out-dir "$OUTPUT_DIR"
```

The provenance JSON contains `entity_id` (or legacy `entity`), integer `version`,
and `sha256`. Conflicting entity aliases are rejected. This
command is documentation only; no final inference has been performed here.

## Limits and next relatedness work

The refit repairs probability support; it does not show a substantial accuracy
gain over the original forest. Distant cousin labels remain difficult and the
class probabilities are not calibrated posteriors. The fixed pedigree templates
and crossover/map stresses cannot certify every possible generating pedigree.
Current input provenance and a final-input run remain separate release steps.

The proposed Task 12 independent-reference refit has now been completed in a
separate artifact directory. It used the six existing calibration replicates
and four new confirmation seeds. Exact accuracy decreased from 37.64% to 37.37%,
so the candidate was not promoted. This corrected the evaluation protocol but
did not improve the method. Whole-family background genealogy, detected-IBD and
consensus validation remain outstanding; see the
[Task 12 outcome](task12_readiness_2026_11_20.md).
