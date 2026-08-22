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

The leading generator-level hypothesis is 20,000 AFR individuals, 30,000 EUR
individuals, a 5,000-generation split, and migration 0.00015, with an
unreported ancestral size near 15,000. This hypothesis accounts for the
roughly 15–20% diversity suppression already measured under the same linked-
selection design in Task 3. Four-arm forward-selection recovery, distance-band
diagnostics, chromosome-arm holdouts, and linked-block uncertainty remain open
gates before any result can be frozen.
