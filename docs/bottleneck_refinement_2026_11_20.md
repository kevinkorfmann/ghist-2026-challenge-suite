# Task 1 bounded coarse-endpoint refinement

Refining every coarse endpoint once with the fine step improves the numerical
evidence without discarding any existing candidate. Across all 53 existing
method problems, the retained union selects a converged best endpoint that
passes the finite-neighborhood check. **52/53 meet every predeclared criterion**;
10 Mb case 20 moments still has only one independent random-start origin near
the best solution. This is a bounded numerical mitigation, not universal
optimizer reliability or biological readiness.

The separate frozen protocol follows the paired comparison described in
`docs/bottleneck_epsilon_extension_2026_11_20.md`. It includes the same 35 spectra
and 53 method problems, without adding newly completed full-length cases or
reading simulation truth. Each of the eight original coarse endpoints per
problem seeds one fine-epsilon refinement: **424 planned refinements**. All 848
original coarse and fine records are copied, hash-verified and retained.
The refinement manifest SHA256 is
`b6c06c91222c09c63b43c179795409b6c9f946fa667b85487bd35b4dc703b3b2`.

Models, bounds, fixed folded likelihood domain, numerical probe source and
per-method production budgets are unchanged from the paired experiment.
Epsilon is 1e-5; moments uses `maxfun=15000, maxiter=500`, dadi uses
`maxfun=500, maxiter=15000`, and both use `factr=1e7, pgtol=1e-5, maxls=20`.
Each initial coarse objective must reproduce within 1e-8 before refinement.
Selection takes the minimum finite objective across original coarse, original
fine and refined endpoints, irrespective of optimizer status. Exact ties retain
an original candidate. Restart agreement uses the best objective for each
original random-start index, so an endpoint and its refinement count once.

The criteria require every planned refinement to complete, no objective
degradation relative to the original coarse/fine union, actual selected
optimizer convergence, at least two distinct original starts within 0.1 of
best, and the same feasible-neighborhood criterion used previously. Retaining
the original union makes objective nondegradation a property of selection;
it does not prove that an optimizer step always improves an endpoint.

Betty array 10851885 completed all 35 cases with one CPU each and at most four
concurrent cases, under a ten-minute limit per case. **424/424 refinements are
finite and complete**, with no exceptions, missing records or timeouts.
400/424 report optimizer convergence; the 24 abnormal endpoints remain in the
record and eligible for objective-based selection. No refinement worsened its
coarse ancestor in this experiment. The summed per-refinement wall time is
537.2 seconds; this excludes setup, original paired fits and Slurm overhead.

| Selected-union diagnostic | Result |
|---|---:|
| Best endpoint reports convergence | 53/53 |
| Best endpoint passes finite-neighborhood check | 53/53 |
| At least two distinct origins within 0.1 of best | 52/53 |
| All predeclared criteria | 52/53 |
| Best endpoint from original fine arm | 27/53 |
| Best endpoint from refinement | 26/53 |
| Best objective improves over original union | 26/53 |

All 17 available full-length method problems meet the numerical criteria;
35/36 shorter problems do so. By method this is 34/35 moments and 18/18 dadi.
The full-length subset still contains only the 11 complete cases available at
the original freeze, not the whole planned biological panel. Every selected
best terminates through relative objective reduction, not projected-gradient
tolerance. The finite probes do not certify stationarity, KKT conditions or
global optimality.

The earlier counterexamples behave as follows:

| Problem | Original coarse/fine-union best objective | Retained refined-union best objective | Distinct origins within 0.1 |
|---|---:|---:|---:|
| 100 Mb case 02, moments | 699.359877595 | 699.359877541 | 8 |
| 10 Mb case 16, moments | 382.013644176 | 376.633413456 | 4 |
| 10 Mb case 17, moments | 434.784550214 | 434.784549899 | 3 |
| 10 Mb case 20, moments | 304.830111360 | 304.830110573 | 1 |

Case 20 remains a failed agreement gate. Origin 2 alone reaches the best region;
the other origins have best objectives approximately 311.298 or 320.745.
Several descendants of origin 2 are close to best, but they cannot establish
independent restart agreement. No further optimizer or threshold search is
performed to clear this failure.

An independent source review confirmed candidate retention, status-independent
selection, exact-tie handling, ancestry counting and actual clipped probe
widths. Synthetic summary checks separately exercise those selection and
ancestry rules. All input/source/output hashes and 424 ancestor links were
verified, and the summary recomputes neighborhood improvements from raw probes.
The reviewer identified a failure-path limit in the frozen diagnostic runner:
if a neighborhood calculation raises after a finite optimizer endpoint, its
single update expression records a failed attempt without that endpoint.
The all-complete gate would then fail. **This path was not exercised: every
refinement and neighborhood completed.** The executed runner is preserved
unchanged; a reusable implementation must preserve finite endpoints before
performing optional diagnostics.

Runtime package versions and source-lock boundaries are those documented for
the paired extension. This is repair development using previously evaluated
spectra, with fixed settings and no truth-based selection. Original scientific
panel outputs remain immutable and are not rescored using these estimates.
The bounded strategy supports a prospective production mitigation with honest
per-arm outcomes, actual selected status and distinct-origin agreement gates;
it does not establish biological performance on new simulations or final data.

Evidence and the reproducible summary are under
`experiments/readiness_2026_11_20/bottleneck/numerical_refinement/`.
Run `.venv/bin/python experiments/readiness_2026_11_20/bottleneck/numerical_refinement/summarize.py`
to audit all planned outcomes and recompute selection and criteria.
