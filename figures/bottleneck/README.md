# Bottleneck figures

Figures are generated from frozen validation artifacts, never from public testing
answers. Every plot is exported as PDF, SVG, and 400 dpi PNG with regular-weight
typography, redundant encodings, and a colorblind-safe palette.

```bash
.venv/bin/python scripts/plot_bottleneck_validation.py \
  --input checkpoints/bottleneck/validation_predictions.npz \
  --output-stem figures/bottleneck/validation
```

The final figure set will include parameter recovery, nuisance-stratified error,
posterior calibration, blockwise stability, component agreement, and observed-data
diagnostics after the real VCF has passed its forensic audit.
