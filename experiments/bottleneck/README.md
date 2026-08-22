# Bottleneck experiment protocol

This screening loop adapts the fixed-budget and auditable-diff principles of
[Karpathy's autoresearch](https://github.com/karpathy/autoresearch) to demographic
inference. It does not modify source code during a run and never reads public
testing answers.

Each candidate receives the same simulation set, grouped validation split,
optimizer seed, epoch count, and batch size. A machine-readable record stores the
complete configuration, source revision, dataset checksum, elapsed time, nominal
and perturbed-domain errors, and 90% interval coverage.

Candidate selection prioritizes normalized nominal-domain RMSE. Smaller penalties
cover overall error, interval miscalibration, and catastrophic performance in any
nuisance stratum. A screen winner is provisional: promotion requires a fresh
simulation shard, five training seeds, ablations, and bootstrap confidence
intervals. Public testing answers remain sealed until the design is frozen.

Promotion also has a chromosome-level gate. Two hundred held-out parameter cells
are each simulated as fifty independent 2 Mb blocks, matching the observed 100 Mb
aggregation exactly. Parameter cells are divided into disjoint calibration,
selection, and final-evaluation sets. An affine latent calibration is fit only on
the calibration cells; the final set is untouched by both fitting and selection.
The inferred-genealogy proxy is diagnostic only until it has a matched translation
validation set, so it cannot displace the summary branch on optimistic tree input.

The fresh promotion used 20,000 new simulations and five seeds. Its grouped
validation score improved from 0.5895 in the architecture screen to 0.5347, with
90% latent coverage of 0.928 and 0.926. The stricter 200-cell chromosome gate did
not justify primary use: on 61 untouched parameter cells, the calibrated summary
branch had median generation relative error 0.234 and median fraction absolute
error 0.030. Its calibrated observed estimate is retained as a diagnostic, while
the final estimator remains mechanistic-dominated pending an independent dadi
cross-check. FastCXT remains ineligible because inferred-genealogy translation
error is not represented by the simulated-tree proxy.

The initial screen varies width, local receptive field, context depth, posterior
mixture size, and regularization. Subsequent screens may vary feature groups,
optimizer settings, genealogy corruption, and nominal-to-robust simulation ratios,
one factor family at a time.
