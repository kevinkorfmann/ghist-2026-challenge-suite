# Relatedness from true IBD segments

Infer meiosis count and relationship labels from exact IBD1 and IBD2 segments.
The current official input contains no HBD segments.

The frozen solution uses sex-map-aware segment spectra and chromosome sharing.
A degree-factorized extension improved balanced accuracy but failed independent
exact and degree confirmation, so the 82.34% exact estimate remains unchanged.
No submission was made.

The November 20 review also provides a separate, checksum-pinned ten-label
[release candidate](../../../configs/task9_current_release.json) that fixes the
historical model's retired-class probability-support failure. Its production
predictor and validator have independent simulated-input checks; see the
[readiness report](../../../docs/task9_readiness_2026_11_20.md). Select it explicitly
with `scripts/predict_task9_current.py --release`. Historical freezes are preserved,
and no final challenge-data inference or submission has been performed.
