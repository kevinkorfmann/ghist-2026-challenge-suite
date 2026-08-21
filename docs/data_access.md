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

