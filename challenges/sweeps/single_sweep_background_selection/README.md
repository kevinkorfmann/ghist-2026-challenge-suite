# Single sweep with background selection

The supported provisional entry point is:

```bash
python -m scripts.predict_task7_sweep_bgs --vcf INPUT.vcf.gz \
  --accessible ACCESSIBLE.bed --dataset testing --raisd /absolute/path/RAiSD \
  --out-dir NEW_OUTPUT_DIRECTORY
```

`blind_callable_median_ten_v1` uses a scan-selected region, callable gap/LD
localization, five composite scales and RAiSD on the same masked SNPs. It retains
raw artifacts and source/input hashes, clips the interval and validates the BED.
It requires chromosome-wide phase; mixed/multiple PS blocks are rejected.

This new method is not scientifically calibrated. The old caller-region,
unmasked archives are preserved but cannot enter the new freeze path. The
freeze CLI requires explicit `--candidate`, `--matched-recovery`, `--bed`, and
`--out` paths and new complete, linked recovery artifacts. No compatible recovery
producer or approved biological release protocol is available. Any future
method package remains `release_ready: false` pending separate scientific review.

See [the readiness review](../../../docs/readiness_sweeps_2026_11_20.md) for the
measured limitations, exact method and next validation dependencies.
