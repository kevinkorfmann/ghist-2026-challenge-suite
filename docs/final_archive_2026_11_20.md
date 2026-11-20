# Final methods and project archive

The November 20 archive covers the work used to prepare all fourteen final GHIST entries: source code, tests, method descriptions, development experiments, selected models, final output files, score records and the submitted methods writeup. The repository remains private.

The [methods schematic](../figures/methods_2026_11_20/ghist_methods.png) follows the visual style of Kevin Korfmann's annual-review figures and gives each task its own method row. The [full submitted writeup](../experiments/final_submission_2026_11_20/writeup_deadline/ghist2026_methods_wiki.md) supplies the detailed model definitions and experimental record.

## Start from the submitted outputs

Install the repository environment using the existing lockfile:

```sh
uv sync --locked --python 3.12 --extra dev --extra references --extra figures
uv run python scripts/verify_final_archive.py
uv run pytest -q
```

The verification command checks the exact SHA-256 of every submitted candidate and matches its scored receipt. It does not fit a model or contact the submission server. The final outputs and receipts are in `experiments/final_submission_2026_11_20/taskN/`. Full numerical scores and valid submission identifiers are in [the final summary](final_scored_summary_2026_11_20.md).

## Code and method map

| Tasks | Procedure | Main code and execution record |
|---|---|---|
| 1–4 | VCF spectra, moments/dadi fits, prescribed model/projection selection and physical scaling | `src/ghist_bottleneck/`, demographic fitting scripts, `experiments/final_submission_2026_11_20/demography/` |
| 5 | Median composite, variant-gap and flank-LD localization | `scripts/predict_single_sweep.py`, sweep method lock and raw prediction JSON |
| 6 | H12 maximum over 50-SNP windows at 12-SNP steps | Soft-sweep prediction/development scripts and `sweeps_takeover/task6_betty/` |
| 7 | Ten-estimator median with multiscale composite, gap, flank LD and RAiSD | `scripts/predict_task7_sweep_bgs.py`, `sweeps_takeover/task7_betty/` |
| 8 | Thresholded RAiSD candidates with spacing and exon snapping | Multiple-sweep predictor, `sweeps_takeover/task8_betty/` |
| 9 | Exact IBD summaries and 600-tree, ten-label random forest | `scripts/predict_task9_current.py`, `configs/task9_current_release.json`, `task9/run_v1/` |
| 10 | Detected-IBD pair classifier and structural pedigree decoding | `scripts/develop_task10_detected_ibd.py`, `experiments/task10_development_2026_11_20/structural_inference_handoff_v4.py` |
| 11 | HBD masking, homolog-sharing and inferred-family context with ExtraTrees | `scripts/develop_task11_pedigree.py`, `scripts/prepare_task11_candidate.py`, corrected context candidate and bindings |
| 12 | Empirical background classifier with the selected top-20 correction | `experiments/submission_package_2026_11_20/task12_deadline/v2/` and its method lock |
| 13–14 | ANGSD GL2/SAF, winsfs spectra and moments demographic fitting | `src/ghist_lowcoverage/`, `scripts/lowcoverage_predict_task13_pragmatic.py`, `scripts/lowcoverage_predict_pragmatic.py`, `source_manifest_v12.json` |

Task 1–4 schema corrections add the official YAML wrapper without altering the estimates. Task 10/12 schema corrections map server-rejected relationship categories to the accepted domain. Original attempts and correction provenance are archived separately.

## Large archives

The GitHub release **ghist-2026-final-20261120** supplies large files separately from the source tree:

- `ghist-2026-work-evidence.tar.gz`: a byte-preserving snapshot of the completed local work and its file manifest, including experiment files too large for the compact Git tree.
- `ghist-2026-frozen-assets.tar.gz`: selected inference source snapshots, trained models, maps, supporting binaries and their bound development assets. Its manifest maps original paths to SHA-256 objects and records 829 checked source/model bindings.
- Companion manifests, checksums and archive receipts allow independent verification.

Download the assets through the release page or authenticated GitHub CLI, then verify checksums before restoring them:

```sh
gh release download ghist-2026-final-20261120 \
  --repo kevinkorfmann/ghist-2026-challenge-suite --dir downloaded-archive
cd downloaded-archive
shasum -a 256 -c SHA256SUMS
cd ..
python scripts/restore_frozen_assets.py \
  downloaded-archive/ghist-2026-frozen-assets.tar.gz runtime --verify-only
python scripts/restore_frozen_assets.py \
  downloaded-archive/ghist-2026-frozen-assets.tar.gz runtime
```

The restoration tool verifies each content object, rejects unsafe paths and refuses to overwrite different existing files. The original runtime root was `/vast/projects/smathi/cohort/korfmann/ghist2026`; the manifest records paths relative to that root. Several historical execution scripts use absolute server paths and deadline guards. For a replay, use the source snapshot identified by the task's method lock, supply the official inputs at its documented paths, and reproduce the recorded runtime environment. Do not treat deadline orchestration wrappers as generic current-date commands or overwrite archived predictions.

Raw organizer VCF/BAM/reference downloads, virtual environments, caches, authentication material and unrelated private messages are outside the publication package. Official input versions, identifiers and hashes remain in the data and task manifests so the inputs can be obtained from their original source. Historical experimental files are snapshots, while the top-level source and tests provide the maintained checkout.

## Verification before publication

The complete test suite passed 1,494 tests. Ruff checks passed for the maintained source, scripts and tests. All fourteen candidate hashes and scored receipts match. The local-work snapshot contains 3,072 files (558,180,506 bytes); large files are included in the release asset and indexed in the snapshot manifest.
