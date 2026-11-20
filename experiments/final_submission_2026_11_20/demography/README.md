# Tasks 1–4 final best-effort handoff — 20 November 2026

The four point predictions in `candidates/` are the definitive outputs of the fixed final-data plan. All 48 planned fit cases reached a terminal state: 34 produced finite, provenance-verified fits and 14 failed numerical integration checks. No cases were missing, censored, or running when the final package was assembled. No external submission was made by this runner.

These are best-effort entries under the user’s deadline instruction. They do not pass scientific release qualification. Convergence, recovery accuracy, misspecification robustness and confidence-interval coverage are separate questions; no intervals were generated. In particular, the selected Task 2 endpoint and two selected Task 4 endpoints report optimizer `ABNORMAL` termination. They remain selected because the locked rule uses the best finite likelihood without hiding lower-objective, unconverged candidates.

| Task | Planned cases | Finite completed | Failed | Selected optimizer convergence |
|---|---:|---:|---:|---|
| 1 | 2 | 2 | 0 | true |
| 2 | 16 | 13 | 3 | false |
| 3 | 6 | 6 | 0 | true, true, true |
| 4 | 24 | 13 | 11 | true, false, false |

Task 1 uses the folded two-epoch moments fit with 16 original starts and the versioned paired coarse/fine/refinement strategy. The eight-start dadi run is retained as a sensitivity analysis. The primary was fixed before inference; a likelihood difference between engines did not change that choice.

Task 2 compares eight separate starts of each existing simple and historical `generator_matched` model on the same unfolded 40×30 spectrum, with the same likelihood mask, five free parameters and optimized theta. The historical model’s fixed hidden epochs remain unverified. Selection is maximum finite multinomial likelihood, with simple-model then origin-order tie breaking. Population 1 is mainland (43 diploids); population 2 is island (24 diploids). The exported migration time is the fitted contact time.

Task 3 uses exponential growth as the fixed primary and the componentwise median of physical parameters from projections 84, 60 and 40. All three primary fits and all three sudden-growth sensitivity fits completed. The alternative convention does not change the predeclared primary. Background selection and convention-dependent identifiability remain scientific limitations.

Task 4 selects the maximum finite likelihood separately within each of the 50×40, 70×50 and 94×68 projections, then takes componentwise physical medians. These are eight original starts per projection. Likelihoods are never compared across projections, and the median is not treated as a single mechanistic model vector. AFR precedes EUR (47 and 34 diploids). No legacy learned calibration was applied.

All 14 failed fit cases failed with the retained reason `Spectrum integration failed after 7 finite refinements: Expected spectrum would exclude observed likelihood bins`. The fixed observed likelihood domain was preserved. Their complete failure records and raw logs are included; failures were not silently dropped from denominators.

Input identities are the authenticated current Synapse versions: Task 1 syn74933318 v3, Task 2 syn74933326 v3, Task 3 syn75193300 v5 and Task 4 syn75193301 v5. Exact MD5, SHA256 and byte counts are in `input_provenance/`. Root separately confirmed current coding BED syn74933457 v1 and recombination maps syn75318370 v2 / syn75318369 v1 against the actual cached bytes in `live_reference_confirmation_20261120.json`; this supplement did not alter the frozen input bindings.

Tasks 1 and 2 use 100,000,000 bp of header exposure. Tasks 3 and 4 use exactly 61,965,066 bp in the fixed noncoding 1–100 kb distance strata, computed from the actual contig lengths and coding BED. The same exposure is passed to each fitter and checked in the returned fit. Tasks 2–4 use valid ancestral-state orientation. Task 1 is folded and has no incomplete genotype records; Task 2’s fully-called fraction is 1.0. Spectrum identities, full preflight audits, remote archive locations and hashes are recorded in `input_spectrum_inventory.json`.

The first Task 1 extraction failed an overstrict unique-position audit. The corrected preflight permits sorted nondecreasing VCF positions and retains every original record. It audits 13,226 adjacent repeated positions in Task 1, with no extra coordinate shift, deduplication, allele change or genotype replacement. The failed job10859290_1, unattempted cancelled dependent array10859293, initial protocol/runner and failure logs are preserved under `failed_extraction_v1/`. Numerical model settings and seeds were unchanged by this harness correction. Successful extraction array10859396 and fit arrays10859397–10859400 are fully accounted for.

The source manifest contains 247 project source/configuration files from source commit0ced7ae72c73199da8fcc49996ddd2d079e67106 (the tested source subsequently merged via PR10). The original protocol bytes/hash are bound across extraction and fitting. Every imported `ghist_` module is checked by its actual resolved `__file__` and SHA256 against the frozen source snapshot, before/after fitting. Input bytes, spectrum hashes and source identities are rechecked before completed results become eligible. All 247 local source hashes and all input bindings were independently rechecked during handoff assembly. Selected raw-fit hashes and final YAML hashes also agree. The packager SHA is 07e075f99be0a58452d95dd4988c6991a2430c4337d9c732f0274971736d25aa.

These checks establish provenance, arithmetic consistency, a fixed selection procedure and output-schema validity. They do not establish that a selected finite endpoint is a global optimum or stationary point, nor that its parameters are biologically accurate. Existing scientific recovery failures remain unchanged; no known answers or held-out challenge truths were used for this final run.

`independent_reviews/` preserves separate saved-value reviews. The Tasks 1/3 and Task 4 reviews refer to earlier immutable packages whose selected candidate bytes are identical to this definitive package. Task 2 is reviewed against this definitive package. The reviews recompute arithmetic and selection from saved outputs only; they do not run additional model evaluations.

Reproduction requires restoring the frozen source files against `source_manifest.json` and the exact raw inputs/spectra against their recorded hashes in a fresh directory. The included Slurm wrappers are historical invocations and contain the original remote working path; they are not relocatable commands. Use a new path-adapted batch wrapper with the same environment and resource settings, preserving the frozen runner/protocol/source/input bytes, and retain the archived outputs unchanged. The runner refuses existing output artifacts. No reproduction job is launched by this archive.

Each job used one CPU and 5,632 MiB, with a 20-minute bound and at most 12 concurrent fit jobs. `terminal_slurm_accounting_utc.psv` was explicitly collected with `TZ=UTC`; the other accounting file preserves the scheduler’s native timestamp presentation. `terminal_audit.json` records the actual UTC collection time and full planned-case denominator. Raw final VCFs and NPZ spectra are intentionally omitted from this compact repository bundle; their immutable locations, byte counts and hashes are retained.
