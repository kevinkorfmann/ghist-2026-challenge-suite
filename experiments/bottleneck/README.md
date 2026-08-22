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

The initial screen varies width, local receptive field, context depth, posterior
mixture size, and regularization. Subsequent screens may vary feature groups,
optimizer settings, genealogy corruption, and nominal-to-robust simulation ratios,
one factor family at a time.
