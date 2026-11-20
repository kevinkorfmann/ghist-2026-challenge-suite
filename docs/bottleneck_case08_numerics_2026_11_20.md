# Task 1 case 08 optimizer diagnostic

The original 100 Mb case 08 selected moments and dadi fits both stopped with
`ABNORMAL` line-search status, despite six converged restarts out of eight for
each solver. This was not an evaluation-budget stop: the selected fits used
159 and 135 evaluations. Their finite estimates and solver agreement did not
justify clearing the convergence gate.

An isolated diagnostic reused exactly that spectrum, folded likelihood domain,
parameter bounds and numerical models. It fixed eight new log-uniform starts
from seed 2026113008 before execution, applying the identical starts to moments
and dadi with forward log-parameter difference steps of 1e-3 and 1e-5 and a
bounded Powell comparison. No simulation truth entered initialization,
comparison or selection. The diagnostic used `maxfun=1000`, `maxiter=500` for
both L-BFGS variants; this differs from production's dadi `maxfun=500` and the
moments wrapper's default `maxfun=15000`; dadi's default `maxiter` is 15000.
Thus its paired epsilon comparison is
controlled, but it is not a repeat of the original optimizer budgets.

Betty job 10851187 completed all 48 planned fits in 107 seconds on one CPU, with
236,788 KB maximum reported RSS. No fit was censored. Original input and source
hashes were verified after execution. The original selected objectives
recomputed exactly; all original panel outputs remain unchanged.

At the original selected estimates, the coarse forward difference has the
wrong sign for the decline-fraction coordinate in log space. Moments gives
+29.923 versus a centered 1e-5 derivative of −54.880; dadi gives +15.020 versus
−69.627. Centered 1e-4 and 1e-5 derivatives agree closely. The large local
curvature makes a 1e-3 forward step substantially biased. This directly supports
finite-difference truncation bias as a cause of the local line-search trouble;
it is not an objective-sign or parameter-transform error.

| Solver / optimizer | Reported convergence | Best negative log likelihood | Worst negative log likelihood | Best centered gradient infinity norm, h=1e-5 |
|---|---:|---:|---:|---:|
| moments, forward 1e-3 | 6/8 | 507.857961 | 507.863676 | 57.008 |
| moments, forward 1e-5 | 8/8 | 507.841506 | 507.841506 | 0.822 |
| moments, Powell | 8/8 | 507.841504 | 21017.320543 | 0.00362 |
| dadi, forward 1e-3 | 7/8 | 508.212862 | 20165.708785 | 30.245 |
| dadi, forward 1e-5 | 7/8 | 508.209035 | 20164.948400 | 0.504 |
| dadi, Powell | 8/8 | 508.209034 | 17905.425444 | 0.00326 |

The selected best fine-step fits report convergence for both solvers, through
relative objective reduction, not the projected-gradient tolerance. Their raw
centered gradients remain about 0.5–0.8, so they are not certified stationary.
Their objectives agree with the best Powell results within 2.084e-6 and
8.869e-7, and their two fitted parameters differ by less than 7.1e-6 relative.
The numerical change is small in likelihood, while materially improving the
gradient approximation and local termination behavior.

Neither a convergence flag nor Powell is a safe automatic fallback. Powell
reports success at several very poor, nonstationary interior solutions. Fine
differences do not remove dadi's poor boundary basins or its remaining abnormal
restart. The diagnostic's boundary-clipped secant slopes are not centered or
KKT tests; its `centered` flag only indicates that both sides existed. The
selected best fits discussed above are interior and use symmetric differences.

This single spectrum was selected because its optimizers terminated abnormally.
It supports a finer-gradient candidate for a broader fixed comparison, not a
global optimizer guarantee or biological readiness claim. The subsequent
paired extension included all 11 complete 100 Mb spectra available at its
freeze plus all 24 earlier independent 10 Mb spectra, using production-matched
budgets and numerical criteria independent of truth errors. It improved 52/53
best objectives but materially worsened one moments problem and failed at least
one criterion in four problems. See
`docs/bottleneck_epsilon_extension_2026_11_20.md` for the full denominator and
counterexamples. A blanket fine-step default change is not supported; no
production default has been changed from these diagnostics.

The executable protocol, all starts and results, independently verified output
hashes, Slurm accounting and reproducible summary are under
`experiments/readiness_2026_11_20/bottleneck/numerical_case08/`.
