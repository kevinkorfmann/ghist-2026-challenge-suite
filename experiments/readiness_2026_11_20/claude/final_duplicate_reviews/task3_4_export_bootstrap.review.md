## Scope and method

Read all 13 listed files once each (13 Read calls, no other files, no commands, no edits). All findings below trace to lines in these files; where a conclusion depends on artifacts not in scope (block archives, `spectrum.py`, audit JSONs), I say so and give the check rather than asserting the outcome. Diploid sizes and chromosome-arm/pedigree independence units assumed as stated.

---

## Finding 1 — The published interval is the sampling distribution of a *different estimator* than the published point estimate

**Severity: High (invalid scientific claim; traceable in code).**

**Where:** `scripts/freeze_growth_bgs.py:16,72-76,109,128`; `scripts/freeze_split_migration_bgs.py:26,65-69,99,120`; `src/ghist_inference/bootstrap.py:259-263`; `scripts/bootstrap_growth_bgs_blocks.py:22,59`; `slurm/growth_bgs_block_bootstrap.sbatch:32`; `scripts/bootstrap_split_migration_bgs_blocks.py:29,53`.

**Cause:** `testing_prediction` is a componentwise **median across three projections** (84/60/40; 50×40/70×50/94×68). Every bootstrap replicate is fit at **one** projection: `require_summary` *requires* `settings['projection'] == 84` (or `[50,40]`) and the producers never fit the other projections. The interval written as `linked_block_interval_95` is therefore the 2.5/97.5 percentile spread of the single-projection estimator, attached to a three-projection median. Projection to 40 vs 84 haplotypes changes both bias and variance (different site counts, different sensitivity to BGS-distorted high-frequency classes), so the median-of-projections estimator has a different sampling distribution — and nothing anywhere checks that the point estimate even lies inside the published interval.

**Verification:** For each parameter, test containment of `testing_prediction[p]` in `linked_block_interval_95[p]` in the frozen JSON. Then take ~20 existing bootstrap replicates, refit the two non-primary projections on the *same* resampled block sets, and compare the percentile width and center of the median-of-three estimator to the primary-only interval. If widths or centers differ materially, the interval must be recomputed per-replicate on the aggregation actually published, or the point estimate must be changed to the primary projection.

---

## Finding 2 — "95" interval propagates only resampling variance, while the freeze gates admit up to 10–15 % systematic recovery error

**Severity: High (numerical convergence vs biological accuracy conflated).**

**Where:** `scripts/freeze_growth_bgs.py:99-105,128`; `scripts/freeze_split_migration_bgs.py:77-88,120`; `src/ghist_inference/bootstrap.py:290-294`.

**Cause:** The freeze passes when known-truth recovery error is below 0.10 (`generations`), 0.05 (`post_growth_scaling`, sizes, split time) and 0.15 (`migration_rate`) — i.e. the method is certified only to ~5–15 % accuracy on simulated truth. The exported interval is purely the block-bootstrap percentile range (`validated_quantile_intervals`), which contains no term for that demonstrated systematic error, nor for BGS-induced model misspecification. The artifact nonetheless publishes the key as `linked_block_interval_95`. The honest flag `interval_coverage_validated: False` is recorded next to it (`freeze_growth_bgs.py:134`), which makes the two statements mutually contradictory rather than qualified: a key asserting a 95 % interval and a field denying coverage validation cannot both stand.

**Verification:** Compute, per parameter, `half_width / point` from the frozen interval and compare with `max |relative_error|` over `grid["results"]` (task 3) and `recovery["relative_error"]` (task 4). Wherever the bootstrap half-width is smaller than the admitted systematic error, the interval provably cannot have 95 % coverage; rename to `linked_block_percentile_range` or convolve the calibrated systematic term.

---

## Finding 3 — `--convention` is a required input whose contents are never read: a dead gate

**Severity: Medium (reproducible code defect).**

**Where:** `scripts/freeze_growth_bgs.py:36,58,64` (and `source_artifacts`, line 158-161); demonstrated by `tests/test_demography_freeze_integrity.py:30`.

**Cause:** `--convention` is `required=True` and loaded into `sources`, but `sources["convention"]` is never dereferenced by any check in `main()`. It only contributes a SHA-256 to `source_artifacts`. The in-scope test fixture passes `'convention': {}` and the freeze succeeds, so an empty JSON object satisfies this evidence slot. The adjacent model-selection gates are also pure trust reads of upstream booleans (`comparison.get("selected_convention")`, `selection_consistent_across_groups`, `smoke.get("preferred")`, `stress.get("freeze_gate_passed")`, lines 84-107) whose producers are not in the hashed source set (`bootstrap.py:25-33`). Convention identity *is* independently enforced for the fits via `require_summary` → `fit['convention'] == settings['convention']` (`bootstrap.py:268-269`), so the defect is a missing gate, not a wrong model label.

**Verification:** Run the freeze with `convention.json` = `{}` and with a file asserting `"selected_convention": "linear"` — both succeed, proving the gate is inert. Then add an explicit assertion on the field(s) that artifact is meant to supply and re-run.

---

## Finding 4 — Bootstrap resampling universe does not equal the fitted exposure denominator, and the two tasks disagree on how blocks are counted

**Severity: Medium-High (denominator accounting).**

**Where:** `scripts/bootstrap_growth_bgs_blocks.py:42-45,50,58`; `scripts/bootstrap_split_migration_bgs_blocks.py:44-46,52-53`; `scripts/freeze_growth_bgs.py:65,152-155`; `scripts/freeze_split_migration_bgs.py:148`; `src/ghist_inference/freeze.py:157-160`.

**Cause:** Three linked mismatches, all in code:

1. Growth drops every terminal partial block (`usable = (ends - starts) == block_size`) and records this (`terminal_partial_blocks_excluded_from_resampling: True`). With 5 Mb blocks that removes up to one block per arm, and the removed material is systematically terminal (distinct recombination rate and BGS regime). The bootstrap variance therefore describes an interior-biased subset, not the dataset the point fit used.
2. Split does the opposite (`within_contig_all_stored_blocks`): it resamples all 1 Mb blocks including partials with equal probability, so unequal-exposure units are mixed and total resampled exposure is itself random. The two published intervals are not the same kind of object.
3. The published exposure denominator is never bound to the denominator used in the physical scaling. `validate_physical_identity` converts θ → N with the fit's own `sequence_length` (`4·µ·L`), but `middle_distance_exposure_bp` is independently reconstructed by summing two audit strata (`noncoding_1kb_10kb + noncoding_10kb_100kb`) while the fits come from the `noncoding_1kb_100kb` stratum. Nothing asserts `sum(strata) == fit['sequence_length']`, nor that the block archive's total exposure matches it. Split has the same unchecked gap via `baseline["physical_exposure_bp"]`.

Secondary: split omits growth's `len(indices) >= 2` guard, so a single-block contig would be silently included contributing zero resampling variance.

**Verification:** From the block archive, print `exposure_bp[usable].sum()`, `exposure_bp.sum()`, and per-contig block counts; compare both sums against each projection fit's `sequence_length` and against the frozen `middle_distance_exposure_bp`. Report the excluded fraction. Any disagreement means the published N's, the published exposure, and the interval rest on three different denominators.

---

## Finding 5 — Split bootstrap is warm-started from the observed point fit, with fewer restarts than the point fit

**Severity: Medium-High (anti-conservative intervals; untested assumption with code mechanism).**

**Where:** `scripts/bootstrap_split_migration_bgs_blocks.py:36,41,51-56` (`initial_scaled=` from `baseline_p50_40.json`), `--restarts 8` (`slurm/split_migration_bgs_block_bootstrap.sbatch:26`) vs growth's 12 (`growth_bgs_block_bootstrap.sbatch:33`); `src/ghist_inference/bootstrap.py:72-73,141-143`.

**Cause:** Every split replicate's optimizer is seeded at the observed baseline optimum. If restart perturbation is local (not verifiable here — `spectrum.py` is out of scope), replicate optima are attracted to the observed point, shrinking the bootstrap spread and producing an under-covering percentile interval; the growth bootstrap has no such anchor, so the two intervals differ in meaning in a second way (cf. Finding 4). The warm start is recorded only implicitly via the `initial_fit` hash and is not surfaced as an interval caveat. Separately, `validate_contract` constrains `restarts` only to be a positive integer, so the bootstrap replicates may legitimately explore less than the point fit did; `validate_fit` enforces only that the fit's `restarts` equals the contract's, not that it equals the observed fit's.

**Verification:** Rerun ~10 replicates (same resampling seeds, so block sets are identical) with cold random starts and with `--restarts 12`; compare the 2.5–97.5 width and the per-replicate optimum spread against the warm-started values. Materially wider cold-start spread quantifies the anchoring bias directly.

---

## Reproducible code defects vs untested biological assumptions

- **Code defects (deterministic, demonstrable in scope):** Finding 3 (dead gate, proven by the in-scope fixture); Finding 1's missing estimator/interval binding and missing containment check; Finding 4's unbound exposure denominators and the growth/split block-universe asymmetry.
- **Untested assumptions (mechanism in code, outcome needs data/runs):** Finding 2's implicit claim that resampling variance dominates BGS misspecification and the 5–15 % calibrated systematic error; Finding 5's assumption that warm starting does not contract the bootstrap distribution; the assumption that 5 Mb / 1 Mb blocks within an arm are exchangeable and long enough that inter-block linkage and BGS-induced long-range correlation are negligible — nothing in these files tests block-length sensitivity.
- **Correctly handled (no finding):** the retained-restart denominator is complete and non-cherry-picked (`bootstrap.py:117-139`: all `restarts` rows present and unique, selected row must be the argmax over valid rows, likelihood must equal the selected row's); θ → physical arithmetic is checked at `rel_tol=1e-12` with moments-consistent units (times as `T·2N`, `m = M/(2N)`, `N = θ/(4µL)`, `freeze.py:157-184`); planned denominators come from the submitted array rather than from present files (sbatch contiguity guards; `--planned-replicates`); any failed or duplicated replicate suppresses *all* intervals; artifacts are never overwritten (`FileExistsError`, atomic temp-file writes).

## Unresolved scientific gates

1. **Source/runtime identity does not pin the executed fitter.** `bootstrap.py:22,25-33` derives `ROOT` from `ghist_inference.__file__` and hashes repo-relative paths, including `src/ghist_{stem}/spectrum.py` — not the `__file__` of the imported module. The sbatch files use `PYTHONPATH=src .venv/bin/python`; if the packages are also installed in the venv, the executed fitter may be a different copy from the hashed one. `runtime_versions` (line 53) is only type-checked at verification (`lines 91-93`) — never compared to the verifying environment — and omits the Python version and BLAS/thread identity, even though `workers=1` is the sole determinism control. The hash set also omits both summarizers and both freeze scripts, i.e. the code that computes the published quantiles and intervals. Gate: assert `module.__file__` equals the hashed path, and compare recorded versions to the live environment.
2. **Monte-Carlo resolution of the tail quantiles.** 128 growth replicates (32 × 4) and 64 split replicates mean the 2.5/97.5 percentiles rest on roughly the 2nd–4th order statistics; endpoint noise is substantial and unquantified. Needs a jackknife-over-replicates or larger plan before the endpoints are quoted.
3. **Boundary-constrained replicates can never enter the interval.** `validate_fit` fails any replicate with `post_growth_scaling <= 1` (`bootstrap.py:115-116`) and, via `continuous_prediction` without `nonnegative`, any replicate with `migration_rate == 0`. Since one failure fails the whole plan, intervals are published *only* conditional on no replicate reaching the boundary — a conditioning the artifact does not state, and which is inconsistent with the interval validator explicitly permitting a zero lower bound for `migration_rate` (`require_summary` line 293; test at `test_demography_freeze_integrity.py:141`).
4. **Independence unit never validated.** Within-arm block resampling is asserted, not tested; no check that arms are the independent unit, no block-length sensitivity analysis, no treatment of the pedigree/haplotype unit relative to the declared `haplotypes` field.
5. **Model-selection and stress evidence is self-attested** by upstream artifacts outside the hashed source set (Finding 3).

## Scope limitations

Bounded to the 13 listed files, read once each. I did not read `spectrum.py`, `validation.py`, `integration.py`, `simulation.py` (`require_physical_scaling`, `SIMULATION_SCALING_SCHEME`), the block `.npz` archives, audit/calibration/grid JSONs, or any other review — so optimizer restart geometry (central to Finding 5), the actual exposure numbers (Finding 4), projection/folding conventions inside the fitters, and the real recovery errors (Finding 2) are unverified here and stated as checks, not conclusions. No commands, fits, edits, or submissions were performed; nothing was tuned against any answer. Severities reflect impact on the validity of published intervals and claims, not on code style.
