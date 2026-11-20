# Reproducible setup

Use Python 3.12 and uv 0.10.8 with the committed lockfile:

```bash
uv sync --locked --python 3.12 --extra dev --extra references --extra figures
uv run --no-sync pytest -q
uv run --no-sync ruff check src scripts tests
```

On Linux, source builds of msprime require a C compiler and GSL development
headers (`build-essential libgsl-dev` on Ubuntu). This also applies to ARM64
hosts where binary wheels may be unavailable. The optional `simulation` extra
adds stdpopsim; SLiM, Ped-sim, hap-IBD, RAiSD, Bonsai, and the GPU-specific
FastCXT stack remain external tools managed by their existing Slurm workflows.
The unit suite does not establish that every external tool is installed.

For Betty, use `source /etc/profile` before calling `sbatch` from a non-login
SSH shell. Standard-memory CPU jobs require 5632M per CPU; large-memory jobs
require 15872M per CPU. B200 GPU partitions use 8192M per CPU and require the
partition-specific CPU/GPU ratio. The repaired job scripts follow those rules.
Heavy inference and simulation must run in scheduled jobs, not on login nodes.

Model files produced under an older scikit-learn version must be loaded in their
original environment unless their compatibility has been independently tested.
The November 2026 lockfile is for new installations; it does not retroactively
certify serialized model compatibility.

### macOS editable-import troubleshooting

If installed dependencies import but `ghist_submission` or another local package
cannot be imported, inspect the virtual environment's `.pth` files with
`ls -lO .venv/lib/python*/site-packages/*.pth`. Python 3.12 skips files marked
with the macOS `hidden` flag. This occurred in the local review environment.
Clear that flag on the startup path files, then verify imports in a new process:

```bash
chflags nohidden .venv/lib/python*/site-packages/*.pth
.venv/bin/python -c 'import ghist_submission, ghist_inference'
```

If those flags return, build the environment outside the synced Documents tree
with `UV_PROJECT_ENVIRONMENT=/absolute/local/cache/ghist-venv uv sync ...`.
The review now uses a local cache environment with `.venv` pointing to it; the
previous environment was preserved for rollback. This is an environment repair;
adding pytest-only import paths would conceal the broken command-line installation.
