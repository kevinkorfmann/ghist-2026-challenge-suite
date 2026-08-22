# Split migration with background selection evidence

The testing VCF contains 47 AFR and 34 EUR diploids, or 94 and 68 haplotypes.
Of 514,421 records, 513,470 are fully called, biallelic, and oriented with the
ancestral allele equal to REF; 951 non-biallelic records are excluded. The
1–100 kb noncoding distance band retains 294,804 variants across 61,965,066 bp.
The final VCF remains sealed, and nothing has been submitted.

`baseline_summary.json` records the three-projection isolation-with-migration
fit from Betty job `7751165`. The 50×40, 70×50, and full 94×68 spectra agree
closely. Their median is 16,192 AFR individuals, 22,949 EUR individuals, a
split 4,388 generations ago, and symmetric migration of 0.000146697 per
generation. The full projection returns 16,366, 23,140, 4,360, and 0.000145969.
These are raw mechanistic estimates rather than an accepted prediction.

An initial 20,000/30,000/5,000/0.00015 candidate was retained as a deliberately
broader validation point rather than promoted. Its four-arm simulation recovers
20,687/29,260/4,804/0.0001589, showing that the estimator is nearly unbiased
and that a speculative 20% size correction would be wrong. The observed fit is
therefore inconsistent with this larger candidate.

Distance-band job `7754282` confirms strong model misspecification in the
selected and sparse extremes. Coding sites return 15,461/22,484 individuals,
4,168 generations, and migration 0.0001606, with an ancestral-size scale of
only 7,564. The 2.44 Mb >100 kb tail contains just 3,446 variants and is highly
unstable, so it is excluded. In contrast, the well-powered 1 bp–1 kb noncoding
band returns 18,010/25,918 individuals, 5,028 generations, and migration
0.0001356. All noncoding sequence returns 16,313/23,358, 4,500, and 0.0001455.
This heterogeneity makes known-truth linked-selection recovery essential; no
single observed distance band is being chosen because it happens to look round.

The near-observed four-arm candidate uses an unreported ancestral size of
13,000 and a submitted truth of 16,000 AFR individuals, 24,000 EUR individuals,
a 4,500-generation split, and migration 0.00015. It recovers
16,305/23,674/4,598/0.0001661. This closely matches the observed projection
median of 16,192/22,949/4,388/0.0001467 and is now the leading generator-level
hypothesis. A separate four-point broad recovery grid is running before any
decision about identity mapping, calibration, or freezing.

The four broad truth points and the two candidate points now form a six-point
held-out calibration gate. Across four-arm aggregates, raw recovery errors stay
within 6.6% for either population size, 3.6% for split time, and 12.4% for
migration. Leave-one-truth-point-out comparison selects the uncorrected identity
mapping with mean competition RRMSE 0.0435 and maximum 0.0684. Median-log-bias,
diagonal-affine, and joint-ridge corrections all worsen mean held-out error.
Accordingly, the observed projection median remains unchanged. A 64-replicate,
chromosome-stratified 1 Mb block bootstrap is the next uncertainty gate.

Chromosome-arm job `7766730` fits each arm separately and then repeats the fit
while withholding each arm. Individual arms are heterogeneous, as expected
under the half-scaled recombination map: split estimates range from 3,770 to
5,272 generations. The more relevant leave-one-arm-out fits are stable enough
to bracket the leading generator: AFR 15,349–16,339, EUR 21,827–23,887, split
4,248–4,687 generations, and migration 0.0001406–0.0001585. This argues for
whole-genome aggregation with chromosome-aware uncertainty rather than choosing
the most convenient arm.

The chromosome-stratified 1 Mb bootstrap contains 108 complete blocks and 64
successful replicates. Its 95% intervals are 14,520–17,037 AFR individuals,
20,937–24,109 EUR individuals, 4,108–4,737 generations, and
0.00013896–0.00016216 migration. The leading 16,000 / 24,000 / 4,500 /
0.00015 candidate lies inside all four intervals.

The leading testing prediction is
16,000 AFR individuals, 24,000 EUR individuals, a 4,500-generation split, and
0.00015 symmetric migration. The aggregation is the median of three stable
1–100 kb noncoding projections followed by generator-level rounding. Identity
mapping wins the six-point held-out calibration, so no correction is applied.
A stricter adaptive numerical integration with one-quarter the original time
step changes every reported parameter by less than 0.0006%. All predeclared
gates therefore pass. The method is frozen in `task4_freeze.json`; nothing has
been submitted.
