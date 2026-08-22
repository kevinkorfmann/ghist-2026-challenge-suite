# Growth with background selection evidence

`observed_audit_summary.json` is the compact, versioned intake record for the
testing VCF, coding BED, and Comeron recombination map. It includes source hashes,
record classes, ancestral orientation, coding and recombination strata, and
100 kb diversity diagnostics. The complete window and frequency-spectrum audit
remains on Betty under `artifacts/task3_growth_bgs/`.

The generic VCF audit ran in job `7745884`, the joint resource audit in job
`7745895`, and the compact summary in job `7745897`. The final VCF remains sealed,
and nothing has been submitted.

`stratum_fit_summary.json` records the first mechanistic growth fits from jobs
`7745903` through `7745905`. These independently used 84, 60, and 40 projected
haplotypes, 64 restarts per spectrum, and the exact physical exposure of each
stratum. The whole-genome growth-ratio median is 2.2407, whereas noncoding
positive-recombination quartiles decline from 2.7516 to 1.6011 as recombination
increases. Over the same quartiles, inferred ancestral size rises from 6,682 to
11,945 while inferred modern size remains much more stable. This is diagnostic
evidence of linked-selection bias, not a frozen prediction. Chromosome-arm
consistency and simulation recovery remain required gates.

The refined middle-distance spectrum combines the 1–10 kb and 10–100 kb
noncoding bands, retaining 61,965,066 physical bases while excluding coding,
the most tightly linked 1 kb flanks, and the sparse 2.44 Mb extreme-distance
tail. Jobs `7745927`–`7745929` estimate, respectively, 5,157/1.995 at 84,
5,052/2.010 at 60, and 4,899/2.033 at 40 haplotypes. Their medians are 5,052
generations and 2.010-fold growth. The near-symmetry around 5,000 generations
and close agreement around two-fold growth make 5,000/2.0 the leading testing
hypothesis, but it remains provisional until linked-block and generator-matched
selection recovery pass.

The 128-replicate, chromosome-aware 5 Mb block bootstrap in job `7745935`
places 5,000 generations inside its 95% interval of 3,983–5,445 and 2.0 inside
its growth-ratio interval of 1.738–2.204. The bootstrap used 20 complete blocks;
terminal partial blocks were excluded from resampling, which shifted its center
relative to the all-base point estimate. A 1 Mb repeat is therefore required to
reduce this edge imbalance while remaining conservative for linkage.

The completed 1 Mb repeat in job `7746623` uses 108 complete blocks. Its 128
replicates are correctly centered: medians are 4,956 generations and 1.985-fold
growth, and means are 4,963 and 1.987. The respective 95% intervals are
3,984–5,913 and 1.851–2.122. This removes the edge-driven center shift while
retaining linked-block uncertainty.

Forward-selection validation is pinned to stdpopsim 0.3.0, its Drosophila
`Gamma_H17` deleterious DFE, the `ComeronCrossoverV2_dm6` map, pyslim 1.0.4,
and SLiM 4.3 commit `f872f3d`. The stack was built and verified in job
`7745934`.

`growth_convention_comparison.json` tests exponential growth against an
instantaneous size jump with equal parameter counts on the same middle-distance
spectra. The sudden model has higher multinomial likelihood at all three
projections, by 14.05, 18.81, and 22.43 log-likelihood units at 40, 60, and 84
haplotypes. It returns approximately 2,825 generations and 1.817-fold growth.
This is not promoted because linked selection can change model preference; the
same comparison is running on known exponential forward simulations.

The matched smoke comparison in job `7749383` resolves the apparent conflict.
Data generated under exponential 5,000-generation, two-fold growth recover
5,371/1.986 with the correct model. When intentionally fit as a sudden jump,
the same data yield 2,879/1.808, closely matching the observed sudden fit near
2,825/1.817. Thus the alternative parameters are the expected misspecified
representation of the leading exponential history. The correct exponential
model also has the higher likelihood on the known-truth simulation, although
the 462 kb smoke provides only a 0.24 log-likelihood separation.

`figures/task3_growth_bgs/task3_validation.{pdf,png}` visualizes the current
distance-band bias, middle-distance projection stability, and all 128 one-Mb
linked-block replicates. It is a validation figure rather than a frozen-result
figure and will be regenerated after the forward-selection recovery grid.
