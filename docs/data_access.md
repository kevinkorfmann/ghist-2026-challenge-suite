# Data access

Challenge files require an authenticated Synapse account. Credentials are never
stored in the repository, copied into logs, or embedded in scheduled-job files.

After setting `SYNAPSE_AUTH_TOKEN` in a private shell session, download one
registered input with:

```bash
.venv/bin/python scripts/download_synapse.py \
  --challenge demography.bottleneck \
  --dataset testing
```

The downloader writes a SHA-256 manifest next to the file. `data/` is excluded
from version control and synchronization; large inputs are transferred to the
cluster separately after checksum verification.


## Current versions (November 20, 2026)

The demography and sweep VCFs were corrected to 1-based coordinates in October.
The archived inventory contains August checksums; matching that inventory does
not establish that an input is current. Task 9's final CSV also has a September
revision. See `configs/synapse_versions_2026_11_20.json` for the metadata snapshot,
including access-denied entries that still require authenticated verification.

The downloader now resolves a version first, downloads exactly that version, and
records the entity ID, version, file handle, size, and SHA-256. New downloads go
to `data/current/<family>/<challenge>/<testing|final>/v<version>/`. Existing raw
inputs, sealed holdouts, and freeze records are preserved. Use `--version` to pin
a known version and `--data-root` to choose a storage location. Pass the new
versioned path explicitly to analysis scripts; historical Slurm scripts still
reference their original August inputs.

```bash
uv run python scripts/download_synapse.py \
  --challenge sweeps.single_sweep --dataset testing --version 2
```

Final sweep inference requires chromosome 9 accessibility, recombination maps,
exons, and its mutation rate (1.29e-8). Testing uses chromosome 7 and 1.44e-8.
Both sets of reference entity IDs are recorded in `configs/challenges.yaml`.
