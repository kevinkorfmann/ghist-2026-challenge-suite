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

Run the offline job checks from the checkout before synchronizing a change:

```bash
.venv/bin/python -m pytest -q tests/test_slurm_resources.py
```

This checks Bash syntax without executing jobs, and checks explicit per-CPU
memory against the repository's documented partition policy. It understands
the header forms used here, including spaced options, `-p`, quoted values, and
inline comments. It does not contact Betty, submit jobs, import job entrypoints,
or read datasets. It does not establish runtime dependencies or live scheduler
configuration; jobs without explicit memory requests still depend on partition
defaults.

Before a planned submission on Betty, initialize the environment and select the
project directory:

```bash
source /etc/profile
cd /vast/projects/smathi/cohort/korfmann/ghist2026
test -d logs && test -w logs
```

All current jobs use relative `logs/` output paths. Their in-script `cd` cannot
relocate files opened by Slurm before the script starts. The synchronization
script creates `logs/`; when using another transfer method, ensure that directory
already exists before submitting. An explicit `sbatch --chdir` also sets the
initial working directory. See the [Slurm sbatch documentation](https://slurm.schedmd.com/sbatch.html).

Review launch-time `SBATCH_*` variables and command-line overrides separately:
they can override checked-in directives. The offline checks are limited to
single-component headers and the options exercised by this repository, rather
than being a complete replacement for `sbatch`. The
[2026-11-20 review](../docs/slurm_preflight_review_2026_11_20.md) records coverage
and remaining prerequisites.
