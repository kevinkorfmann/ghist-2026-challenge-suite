# Tasks 2–4: numerical recovery and independent validation

The nine formerly censored neutral runs and the unfinished unscaled split/BGS
fit now finish. The cause was a reproducible non-progressing numerical retry,
not evidence that those parameter configurations required days of computation.
Execution is repaired, but the completed results expose estimation errors and
optimizer failures that still prevent scientific promotion. No final challenge
data were fitted and no predictions were submitted.

This is a new evaluation of the repairs after `b847c4c`. It preserves the
[earlier scientific evaluation](scientific_demography_2026_11_20.md) and its
incomplete outcomes. Public testing answers and `benchmarks/holdout` were not
read or used for this work.

## Numerical failure and repair

Installed moments 1.7.1's adaptive migration integrator restores the previous
spectrum after a negative proposal. Once its timestep floor is reached, some
proposals repeat with the same integration time and timestep indefinitely.
An independent Python trace, without changing the third-party source, observes
16 identical rejected steps for both a Task 2 and a Task 4 parameter proposal.
The [reproducer and exact states](../experiments/readiness_2026_11_20/demography/native_nonprogress.json)
identify the moments source SHA-256. The small trace terminates deliberately;
it does not infer non-progress merely from elapsed wall time.

The wrappers now use a finite sequence of whole-model, non-adaptive
integrations at successively halved timesteps. A retry requires validity in
every likelihood bin and agreement of two adjacent valid spectra. The default
is seven attempts, with relative L1 tolerance 1e-5. Invalid spectra are never
clipped and observed bins are never removed. Failure is explicit if refinement
does not resolve the problem. For the pathological Task 4 proposal, two
negative active bins persist through all seven resolutions, at approximately
−5.08e-8. Rejecting this proposal is appropriate; timestep refinement alone does
not repair every spatial approximation error.

At the fixed fitted parameters, a separate check halved the integration
timestep for all **53 Task 2/4 fitted models**, including the new 10 Mb panel
and the unscaled BGS fit. Maximum absolute log-likelihood change was
**0.015813**, and the largest relative change in a positive observed bin was
**0.1088%**. Recomputed original scores agreed with recorded scores within
2.4e-11. These are timestep sensitivity checks, not global optimizer or
biological accuracy guarantees. The complete per-fit checks are in
[likelihood_refinement.json](../experiments/readiness_2026_11_20/demography/likelihood_refinement.json).

All three fitters now retain the actual optimizer termination message, warning
flag, iterations, function evaluations, selected-restart convergence, and
converged-restart count. Exact parameter-bound starts are accepted for
bootstrap/refit initialization while boundary warnings remain intact. Freeze
paths require reported convergence for Task 2's standard and refined fits,
Task 3's contributing projection fits, and Task 4's contributing projections
and refined fit. Legacy artifacts lacking these fields require regeneration;
a finite likelihood alone cannot pass. Growth summaries now preserve these
diagnostics rather than discarding them. Convergence remains separate from
scientific adequacy and the pre-existing model, provenance and calibration
gates.

An independent cross-review reproduced a summary/gate bypass: two projection-40
inputs, an unconverged estimate of 1,000 generations and a converged estimate of
3,000, yielded a 2,000-generation median but retained only the latter fit's
convergence metadata. Summaries now reject duplicate projection inputs. Freezes
require the complete prespecified projection sets: 84/60/40 for Task 3 and
50×40/70×50/94×68 for Task 4. Partial successful runs cannot stand in for the
planned aggregate.

The same cross-review reproduced a Task 2 freeze accepting a NaN calibrated
prediction despite converged raw fits. Final, standard and refined predictions
now require finite physically valid values, including contact time no earlier
than the split; nominal stress risks must be finite and nonnegative. Valid zero
migration and contact times at zero or at the split are preserved. Strict atomic
JSON writing rejects nonfinite values elsewhere in the frozen record as well.
These export guards change no fitted estimate or scientific evaluation result.

## Completed neutral panels

Every generating unit is a full independently simulated linked chromosome.
Simulation and optimizer seeds, projections, truth points and misspecification
conditions are recorded per replicate. No windows are treated as independent
replicates. Global msprime ploidy is 2. Parameter points and estimator settings
were fixed before examining this evaluation; no model was tuned to these
outcomes.

| Evaluation | Planned / completed | Censored or failed datasets | Purpose |
|---|---:|---:|---|
| Original 2 Mb panel replay | 22 / 22 | 0 | Reproduce all configurations, including seven censored Task 2 and two censored Task 4 runs |
| Task 3 diagnostic replay | 8 / 8 | 0 | Recover optimizer metadata on the same saved spectra; these are not new independent simulations |
| New 10 Mb neutral panel | 28 / 28 | 0 | Four independent replicates per Task 3/4 point and condition, using disjoint new seed streams |

All 13 previously completed 2 Mb records retain **exactly identical parameter
estimates and likelihoods**. All eight Task 3 diagnostic replays likewise match
their original estimates and scores exactly. Completion here means that both
model fits returned and were checkpointed, not that their accuracy or optimizer
termination passed. Four selected fits in the original Task 2/4 replay and
seven in the new panel report L-BFGS-B `warnflag=2`, `ABNORMAL` termination.
These are retained and cannot pass the new freeze guard. Increasing the
iteration limit alone does not address this termination condition.

### Task 2: hidden history remains consequential

The complete 2 Mb results retain the simple secondary-contact fit and the
historical hidden-epoch convention for every replicate. The latter is still
named `generator_matched` in the API; that name does not establish a match to
the current challenge generator. Under the simple matched model, maximum
contact-time error is 65.25% at point 0 and **247.13%** at point 1, despite
optimizer convergence. Under asymmetric migration, the simple model's maximum
contact-time error is 157.88%.

The independently imposed hidden size changes do not match the historical
fixed-epoch assumptions. In that misspecification condition, one historical
model fit estimates migration about **46.6 times** its true value, and the
maximum contact-time error is 341.64%. Boundary/restart diagnostics flag these
hidden-history fits as unsupported. Both historical-model fits to the
asymmetric condition have abnormal selected-restart termination. These results
establish failure strata; they do not justify selecting a favorable model from
the simulation truths. There is no new target-length Task 2 or interval-coverage
validation in this pass.

### Task 3: exposure helps strong growth, but not all weak-growth ambiguity

The following are maximum absolute errors across four independent 10 Mb
replicates for the generating convention, at projection 30. Matched fits are
used to diagnose model recovery; they are not an oracle convention selector
for challenge data.

| Generating condition | Maximum timing error | Maximum growth-ratio error | Matching convention has higher likelihood |
|---|---:|---:|---:|
| Point 0, exponential | 9.93% | 4.45% | 3 / 4 |
| Point 0, sudden | 3.78% | 4.79% | 4 / 4 |
| Point 1, weak exponential | **71.58%** | **53.36%** | 2 / 4 |
| Point 1, weak sudden | 22.42% | 6.68% | 3 / 4 |

Both strong matched conditions and weak exponential have converged selected
restarts in every replicate. One weak sudden matched fit and one alternative
sudden fit to strong exponential data have abnormal termination. All likelihood
comparisons, including these cases, are retained with convergence indicators.
Weak exponential convention differences range from −0.068 to +0.525 log units.
More neutral sequence does not remove the weak-growth identifiability problem,
and these neutral runs do not validate BGS correction or the growth convention
on current challenge data.

### Task 4: projection and misspecification still matter

For the two matched 10 Mb parameter points, maximum split-time errors at
projection 20×16 are **9.67%** and **11.84%**, respectively. At projection 12×10,
they are 9.34% and **24.40%**. At point 1, maximum AFR-size error increases from
18.57% at the larger projection to 30.03% at the smaller projection. The
asymmetric condition has maximum AFR-size errors of 29.86% and 33.64% at the
two projections. The asymmetry changes generating directional migration rates;
the scalar migration comparison uses the recorded nominal rate and is not a
claim that the symmetric fitted model is correctly specified.

Five selected Task 4 fits in this panel have abnormal termination: two at the
larger projection and three at the smaller projection. Restart agreement alone
had labeled all of these numerically supported, demonstrating why actual
termination status must be a separate gate. All estimates and flags remain in
the aggregate; nonconverged cases were not discarded to improve reported error.

## Unscaled BGS split diagnostic

Betty Slurm job **10850252** resumed the exact earlier Q1 tree sequence, verified
against SHA-256
`30c740212a73e15cf8e36eef44470bfa8890198e6192f7805791f9db6f447eb4`.
The simulation was not rerun. Extraction and inference completed in **1m52s**,
exit 0, with moments 1.6.0. Original generation used stdpopsim 0.3.0, SLiM 4.3,
physical sizes and rates, Q=1, burn-in 2, a 1 Mb 2L segment and seed 202611221.
The exact resume command is in
[resume_split_q1.sbatch](../experiments/readiness_2026_11_20/demography/resume_split_q1.sbatch).

The selected fit converged. Of 12 restarts, nine returned finite valid fits,
eight converged, and three were explicitly rejected after persistent invalid
spectra. Signed errors were **−21.44% AFR size, −16.36% EUR size, +47.14% split
time and −22.66% migration**. Fixing the stall therefore reveals substantial
error even in this unscaled diagnostic. One short linked segment with its
particular burn-in cannot certify BGS calibration, convergence of burn-in,
target-scale recovery or uncertainty coverage. Previous independent-seed
Q5/Q10 comparisons remain descriptive; this run is not a paired attribution
of their differences to scaling.

## Resumption, seeds and provenance

Simulation tree checkpoints bind the generating configuration to the tree
SHA-256. Importing an existing tree requires the caller's original generation
arguments and an explicitly recorded tree hash; a hash alone does not prove
those biological arguments. Spectra are saved before inference and bind
metadata to spectrum contents. Completed model fits are checkpointed
individually and reused only with matching spectrum digest, protocol and
settings. Interrupted second fits preserve the first fit and simulation.
Recorded failed attempts can be retried with their failure history retained.
Unknown checkpoint statuses are rejected before calling the fitter or writing,
so a corrupt status cannot silently trigger recomputation and overwrite its
record. Writes are atomic. These checkpoints prevent repeating completed expensive
simulation; they do not resume an optimizer halfway through a restart.

The previous arithmetic seed allocation repeats between configuration 1,
replicate 0 and configuration 0, replicate 100. `counter_v2` allocates separate
ancestry, mutation and two fit streams using a fixed 16-configuration radix.
It preserves existing allocations when adding replicates and rejects overflow.
`legacy_v1` remains only for exact old-panel replay below 100 replicates. The
new 10 Mb panel uses `counter_v2`; its streams are also disjoint from the
specific original-panel streams. Regression tests check more than 100
replicates. Seed allocation is not a claim of statistically independent
hardware randomness; these are distinct deterministic pseudorandom streams.

The [summary](../experiments/readiness_2026_11_20/demography/summary.json) contains
all replicate-wise errors, termination flags, denominators, likelihood
comparisons and old/new equivalence checks. Raw new-panel records retain an
obsolete human-readable `scope` string saying 2 Mb; their authoritative
`length`, simulation specification and saved spectra are **10,000,000 bp**.
The script now formats scope from requested length; the raw run records were
preserved. Source snapshots, runtime versions, exact commands, tests and
artifact hashes are recorded in
[verification.json](../experiments/readiness_2026_11_20/demography/verification.json).
The final source adds input validation, failure-resume metadata and export
guards after some panels were launched; it does not retrospectively claim
those older processes ran the final file hashes.

## Remaining gates

This pass repairs numerical non-progress and makes failures reviewable. Before
promotion, resolve abnormal optimizer termination without changing the
likelihood domain; validate parameter recovery across prespecified prior and
boundary strata at appropriate sequence lengths; test model convention and
hidden-history robustness; regenerate physical-unit BGS calibration with
multiple independent Q1 references and burn-in checks; and evaluate linked
block intervals for empirical coverage. Current-data identity, full extraction,
final predictions and submission validation remain separate deferred work.
All **128 focused tests pass**, including the native-loop regression,
checkpoint interruption/failure recovery, seed uniqueness, physical-unit
contracts, finite-but-unconverged freeze rejection, duplicate/incomplete
projection rejection, the real Task2 freeze entrypoint's nonfinite-input
rejection and unknown-status checkpoint preservation. Critical Ruff checks
pass. Neither these tests nor the completed small panels establish perfect
inference.
