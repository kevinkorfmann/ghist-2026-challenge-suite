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
