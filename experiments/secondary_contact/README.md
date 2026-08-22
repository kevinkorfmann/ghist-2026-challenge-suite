# Secondary Contact evidence

This directory contains compact, versioned evidence for the Task 2 method
selection. Raw VCFs, simulated chromosomes, inferred genealogies, and fitted
weights remain on Betty and are excluded from Git.

- `resolution_p40_p50_32.json`: decisive matched-chromosome comparison between
  the two final projection candidates.
- `resolution_broad_8.json`: matched comparison spanning 40-by-30 through the
  full sample spectrum.
- `resolution_varying_demography_p40_p50_48.json`: independent projection
  comparison across 24 demographic settings with paired 100 Mb replicates.
- `dual_anchor_transfer.json`: point-grouped calibration, selection, and
  untouched evaluation with two independent center anchors.
- `finite_site_bootstrap.json`: exact finite-site bias and conservative
  sampling intervals from 32 chromosomes.
- `stress_comparison.json`: raw recovery under nominal, recombination,
  mutation-rate, local-size, asymmetric-migration, and polarization conditions.
- `stress_transfer.json`: transfer of the frozen correction to those six
  conditions, including every paired replicate score.
- `task2_freeze.json`: accepted testing prediction, hard-gate evidence, source
  hashes, and Betty job provenance.

No result represented here has been submitted.
