# SLURM workflows

All substantial simulations and model fitting run as scheduled jobs on `betty`.
Login-node use is limited to synchronization, environment inspection, and job
submission/status checks.

Remote project root:

`/vast/projects/smathi/cohort/korfmann/ghist2026`

Runtime gates are separated from model fitting. `bootstrap_cxt.sbatch` verifies
the pinned autoregressive translation checkpoint on CPU, while
`bootstrap_fastcxt.sbatch` compiles the pinned state-space kernels and requires a
successful B200 forward/backward pass. Neither runtime gate performs challenge
inference or submission.
