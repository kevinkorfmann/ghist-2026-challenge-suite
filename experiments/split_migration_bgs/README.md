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
