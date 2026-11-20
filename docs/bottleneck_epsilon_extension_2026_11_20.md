# Task 1 paired finite-difference extension

The fixed paired comparison does **not** support changing both production
optimizers directly from a forward log-parameter step of 1e-3 to 1e-5.
The finer step improves 52 of 53 best objectives and most local numerical
diagnostics, but materially worsens one moments result. Four comparisons fail
at least one predeclared criterion. Production defaults and the original
scientific panel outputs remain unchanged.

The manifest froze every complete full-length spectrum then available: cases
00–10 of the planned 24-case 100 Mb panel. Cases 11–23 were unavailable and were
not added later. All 24 earlier independent 10 Mb spectra were also included.
This yields 35 spectra, 53 method problems (moments throughout and dadi on the
matched scenarios), eight fresh random starts per problem and two epsilon
arms: **848 planned fits**. The full-length subset is conditional on availability
and may favor computationally easier configurations. This is a same-spectrum
numerical development comparison, not new biological replication or a new
recovery/confirmation panel.

The frozen manifest SHA256 is
`8c2e403832193cb437e472baecc371aa9f4e6809122df5a027b0d136c16ae7ba`.
Seed root 2026113010 and `SeedSequence([root, dataset_code, case_index])` generate
eight log-uniform starts, identical between arms. No simulation truth enters
initialization, optimization, comparison or selection. Both arms use the same
folded SFS, fixed likelihood domain, bounds, models and production-matched
optimizer controls. Moments uses `maxfun=15000, maxiter=500`; dadi uses
`maxfun=500, maxiter=15000`. Both use `factr=1e7, pgtol=1e-5, maxls=20`.
Dadi grids are sample size plus 20, 30 and 40. The only optimizer-arm difference
is epsilon. A superseded preflight manifest had dadi's maxiter set to 500;
it was corrected and archived before any execution, with unchanged inputs and
starts. No runs used that preflight manifest.

Each problem must retain all planned outcomes, avoid loss of finite fits,
avoid a best-objective deterioration greater than 1e-3, select a best finite
fit that reports optimizer convergence, have at least two starts within 0.1
log likelihood of that best, and have no tested feasible neighbor improve
negative log likelihood by more than 1e-4. Selection is strictly by best finite
objective, including unconverged candidates. The neighborhood tests all eight
nonzero directions in `{-1,0,1}²` at log steps 1e-3, 1e-4 and 1e-5, clipped to
the parameter bounds. Actual step widths are recorded; only equal nonzero
widths are labeled centered. These finite probes do not certify KKT conditions,
stationarity or a global optimum.

Betty array 10851420 completed all 35 jobs, each using one CPU, with at most
four concurrent jobs and a ten-minute allocation. Individual elapsed times
were 16–176 seconds. **848/848 fits are finite and complete**, with no missing,
exception or censored records; unconverged optimizer terminations remain
explicitly recorded. All 35 completion records and output hashes were verified.

| Numerical measure | Coarse 1e-3 | Fine 1e-5 |
|---|---:|---:|
| Reported converged fits | 377/424 | 405/424 |
| Selected best reports convergence | 47/53 | 52/53 |
| Selected best passes finite-neighborhood check | 11/53 | 52/53 |
| At least two independent starts within 0.1 of best | 51/53 | 51/53 |

Fine best objectives improve in 52/53 problems. Fine-minus-coarse differences
range from −0.044605 to +14.891253, with median −0.001465. All predeclared
criteria hold for 49/53 problems: 16/17 available full-length problems and
33/36 shorter problems. By method, this is 31/35 moments and 18/18 dadi
problems. Every converged selected best terminated through relative objective
reduction; none of the selected best fits terminated through projected-gradient
tolerance. Better convergence counts cannot override the following failures.

| Counterexample | Coarse best objective | Fine best objective | Failed criterion |
|---|---:|---:|---|
| 100 Mb case 02, moments | 699.360396 | 699.359878 | Selected fine best stops `ABNORMAL` |
| 10 Mb case 16, moments | 382.013644 | 396.904897 | Objective worsens by 14.891253; neighbor improves by 0.028670 |
| 10 Mb case 17, moments | 434.788293 | 434.784550 | Only one random start within 0.1 of best |
| 10 Mb case 20, moments | 304.839318 | 304.830111 | Only one random start within 0.1 of best |

The case 16 regression is substantial: coarse restarts 1, 4 and 7 reach the
better region, while their paired fine runs stop in worse regions. All eight
fine fits report convergence, including the selected best, which still has an
improving feasible neighbor. This demonstrates why a convergence flag alone
is insufficient and why the favorable aggregate cannot justify a blanket
epsilon switch. The selected fine full-length case 02 has eight starts within
0.1 and no improving tested neighbor, but its actual abnormal status remains
a failed convergence gate. Cases 17 and 20 improve locally yet retain poor
restart agreement under both arms.

An independent agent checked the frozen source, manifest, 35 input hashes,
objective transformation, paired controls and boundary-probe definitions.
The summary also recomputes neighborhood improvements and flags from all raw
probe values, checks optimizer status/flag consistency, verifies output hashes,
and accounts for every planned start. Separate quadratic boundary checks
confirm that asymmetrically clipped probes are not called centered and that
a feasible bound minimum passes without being described as stationary.

Runtime versions are moments-popgen 1.6.0, dadi 2.4.4, SciPy 1.16.3 and NumPy
2.4.6. Runner sources, inputs and the shared likelihood-domain guard are
hash-locked. Runtime package versions are recorded, but package implementation
files and the imported atomic-checkpoint helper are not source-hash-locked by
this manifest. The standalone objective is equivalent to the production model
and domain at the evaluated endpoints; the experiment is not an exact replay
of the original production optimizer trajectories.

All evidence and the reproducible summary are in
`experiments/readiness_2026_11_20/bottleneck/numerical_extension/`.
Recompute with
`.venv/bin/python experiments/readiness_2026_11_20/bottleneck/numerical_extension/summarize.py`.
The original narrow case 08 evidence is separately retained under
`numerical_case08/`. A subsequent, separately frozen mitigation diagnostic may
refine each coarse endpoint using the finer step, retaining every original
candidate and counting agreement by original random-start ancestry. It must
be evaluated before any production change; it would be repair development,
not independent biological confirmation.
