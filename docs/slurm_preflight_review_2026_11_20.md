# Slurm preflight review — 2026-11-20

This read-only review covered the 149 checked-in `slurm/*.sbatch` files. It found
no demonstrated Bash syntax error, missing literal repository entrypoint, or
explicit per-CPU memory mismatch. No job was submitted or executed, no remote
state was changed, and no dataset or result artifact was read. Generated job
scripts under `experiments/` were outside this review.

| Check | Result |
| --- | --- |
| Bash syntax (`bash -n`) | 149/149 pass |
| Explicit per-CPU memory | 133/133 match the repository's partition policy |
| No explicit memory option | 16 jobs rely on scheduler defaults |
| Literal `python scripts/*.py` invocations | 144 invocations; all 96 distinct files exist |
| Literal `python -m ghist_*` invocations | 14 invocations; all 5 distinct module sources exist |
| Scheduler output paths | All 149 use relative `logs/`; none specifies `--chdir` |

The actionable defect was in the existing resource regression test: it only
recognized `--partition=...` and `--mem-per-cpu=...`. Spaced options and `-p`
could bypass the check, and directives after the first executable line were
incorrectly considered active. The revised tests cover those cases, preserve
quoted hashes and inline comments, normalize equivalent memory units, reject
conflicting memory modes, and run non-executing Bash syntax checks with a clean
startup environment. This improves regression coverage; it does not change any
job's scheduler behavior.

Slurm processes only the initial directive block, treats shell variables there
literally, and lets environment/command-line options override script settings.
Header semantics were checked against [SchedMD's sbatch documentation](https://slurm.schedmd.com/sbatch.html)
and its [header parser](https://github.com/SchedMD/slurm/blob/master/src/sbatch/opt.c).
Memory normalization follows the integer units and KiB rounding in
[SchedMD's argument parser](https://github.com/SchedMD/slurm/blob/master/src/common/proc_args.c).

Two launch prerequisites remain. Submit with the project root as the initial
working directory and a pre-existing writable `logs/` directory; a subsequent
shell `cd` cannot repair initial log-file resolution. The existing synchronization
script creates that directory. Also verify the effective partition memory and
CPU/GPU policy on Betty for jobs that inherit defaults; omission itself is valid
and was not classified as a failed launch. The jobs relying on memory defaults
are:

- `genoa-std-mem`: `bottleneck_chromosome_evaluate`,
  `bottleneck_chromosome_validation`, `bottleneck_fastcxt_dataset`,
  `bottleneck_merge_search`, `bottleneck_moments_chromosome_evaluate`,
  `bottleneck_sfs_bootstrap`.
- `b200-mig45`: `bootstrap_fastcxt`, `bottleneck_experiment_smoke`,
  `bottleneck_experiments`, `bottleneck_fastcxt_observed`,
  `bottleneck_fastcxt_train`, `bottleneck_mamba_smoke`, `bottleneck_promotion`,
  `bottleneck_select_screen`, `gpu_smoke`.
- `dgx-b200`: `bottleneck_train`.

Names above omit the `.sbatch` suffix. There is no proposed job-script fix without
evidence that the live defaults or launch directory are wrong. The entrypoint
inventory establishes source existence only: dynamic shell targets, installed
modules and tools, input files, CLI semantics, dependencies, runtime resources,
and scheduler acceptance remain unverified. The resource helper is scoped to
single-component headers and the option forms used here, not the complete Slurm
CLI grammar.

Validation: `pytest -q tests/test_slurm_resources.py` passed 318 checks, including
149 resource-policy cases and 149 Bash syntax cases. `ruff check
tests/test_slurm_resources.py` passed.
