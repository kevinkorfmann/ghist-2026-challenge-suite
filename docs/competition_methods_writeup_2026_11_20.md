# GHIST 2026 methods — ancestrum

## Author

Kevin Korfmann (Synapse: ancestrum).

This writeup describes the fixed methods used for the fourteen challenge entries. Most entries were submitted as best-effort methods despite unmet internal scientific qualification criteria. Successful execution, valid formatting and complete output coverage are distinct from biological accuracy. Methods were selected using recorded simulation/development comparisons and fixed inference policies; final scoring feedback was not used to refit these models. Rejected submissions received documented schema corrections while their original files were preserved.

## Task 1 — Bottleneck

A folded site-frequency spectrum (SFS) was extracted from the supplied VCF and fitted to a two-epoch population-size model using moments. The mutation rate was 3.42×10⁻⁸, with callable exposure taken from the input contract and contig lengths. The fixed primary optimizer policy paired coarse and fine numerical steps with refinement over sixteen starts; dadi, with eight starts, supplied an independent sensitivity fit. The moments estimate was primary, with dadi available only if moments failed to produce an eligible result. Scaled parameters were converted to generations and post-decline population fraction using the fitted mutation parameter. Numerical sensitivity had required substantial development work; this entry does not establish universal parameter recovery or interval coverage.

## Task 2 — Secondary contact

The ancestral-allele annotations and mainland/island sample labels defined an unfolded joint SFS, projected to 40×30 chromosomes. Globally monomorphic corners were excluded while fixed differences between populations were retained. moments fitted both a simple isolation-then-contact model and a historical model with specified hidden size epochs, using mutation rate 8.4×10⁻⁸. The prespecified choice was the highest finite multinomial likelihood on the same SFS among the two five-parameter models and eight starts per model; theta was optimized in both. Invalid model spectra were rejected rather than clipped. Historical hidden epochs remain assumptions, and the selected fit reported abnormal optimizer termination. The finite point estimate is therefore an explicitly uncertain best-effort result.

## Task 3 — Growth with background selection

We used the unfolded SFS from noncoding sequence 1–100 kb from coding intervals, with its corresponding physical exposure, to reduce direct coding-site effects. moments fitted exponential growth at projections 84, 60 and 40, with eight starts per projection and mutation rate 3.36×10⁻⁸. The reported physical parameters were the componentwise medians of all three fixed exponential fits; sudden growth was a sensitivity analysis rather than an outcome-selected replacement. These filters do not remove background selection. Independent simulations showed weak-growth identifiability and selection/model-convention sensitivity, so no calibrated background-selection correction or confidence interval is claimed.

## Task 4 — Split and migration with background selection

AFR/EUR samples and ancestral alleles defined joint spectra from the same 1–100 kb noncoding distance strata, using exact stratum exposure. moments fitted isolation with migration at projections 50×40, 70×50 and 94×68, with mutation rate 1.015×10⁻⁸. Each projection used the best finite likelihood from eight fixed starts; the final physical parameters were componentwise medians across the complete three-projection set. Integration used bounded timestep refinement and fixed likelihood support. Two of the three selected projection fits did not report convergence. Projection sensitivity, migration asymmetry and background selection remain material limitations; raw likelihoods from different projections were not directly ranked.

## Task 5 — Single sweep

The fixed three-estimator callable-sequence method combined composite-scan, gap and linkage-based localization, with the search region determined by the scan rather than the simulated selected position. The median position received a 60 kb half-width. Final coordinates were moved or clipped only according to the supplied callable mask and chromosome bounds. The final estimator span exceeded the earlier 50 kb agreement criterion, so this is not a consensus-qualified prediction. Earlier short-region simulation recovery did not establish chromosome-wide coverage. RAiSD was not invoked by this task's final method.

## Task 6 — Soft sweep

The submitted method used the previously development-selected H12 haplotype-homozygosity scan with 50-SNP windows and no minor-allele-frequency cutoff. Its selected position received a 75 kb half-width, subject to the common callable-mask and chromosome-boundary rules. This was a fixed development choice without fresh independent confirmation. Weak or incomplete sweeps, alternative soft-sweep mechanisms and whole-chromosome localization remain important failure domains; no general coverage guarantee is attached to the interval.

## Task 7 — Single sweep with background selection

The fixed masked ten-estimator localization method included RAiSD with 50-SNP windows and used its prescribed median position, with an 80 kb half-width. The callable mask governed admissible positions and interval clipping. The exact ten-estimator workflow lacked sufficient matched biological recovery, and a 500 kb composite estimator disagreed with the final aggregate. Those limitations were retained rather than used to select a different final estimate. An initial RAiSD runtime-library failure was repaired without changing scientific settings or inputs.

## Task 8 — Multiple sweeps with background selection

RAiSD peaks in fully callable windows were screened using a provisional threshold: twice the largest callable Task 7 score more than 2 Mb from that task's highest peak. Eligible peaks were processed in descending score order, with lower-position tie breaking and at least 1 Mb separation. Peaks could be snapped to callable exon bases within 60 kb; 60 kb half-width intervals were clipped and overlapping intervals merged. The frozen fallback was one highest callable peak if no background reference or eligible candidate existed. The actual run used the completed Task 7 reference and produced two intervals. This transferred background is not a calibrated null, and the method does not establish accurate sweep-count inference. The exon annotation was interpreted as one-based inclusive and converted to BED; that convention was not explicitly stated by the challenge.

## Task 9 — Relatedness from IBD segments

Exact pairwise IBD1/IBD2 segments were summarized by 218 fixed features, including lengths, counts, distributional and chromosome summaries, IBD2 proportions and sex-specific map-derived lengths. A scikit-learn random forest used 600 trees, balanced-subsample weights, maximum feature fraction 0.8 and minimum leaf size two, trained on 60,266 pairs from six Ped-sim calibration replicates. A fixed ten-label release was tested on two subsequent independent simulation datasets, with 89.10% and 87.20% exact-label accuracy. These figures are simulation accuracy, not the challenge's weighted Jaccard score. Final inference covered all 997 unordered pairs present in the official segment file. No absent zero-IBD pair was invented; distant-cousin discrimination and uncalibrated class probabilities remain limitations.

## Task 10 — Relatedness in a large human family

Phased genotypes were restricted to the frozen human map domain and processed with hap-IBD. Detected haplotype edges were converted to diploid sharing states, requiring two distinct homolog matches for IBD2. A frozen pair classifier supplied relationship probabilities; a neighbor-joining structure based on probability-derived meiosis distances was refined by bounded positive-integer edge moves and scored over admissible roots using composite pair probabilities. A separately frozen whole-cohort pair-only fallback applied only when the initial tree had no accepted root; it was not needed on the final input. All 34,453 pairs of 263 samples were emitted. Fresh large synthetic cohorts gave 67.77% pooled exact-label accuracy, below the prospective 70% practical gate. Following an unscored schema rejection, 3,040 unsupported first-cousin-thrice-removed labels were merged into the server-accepted second-cousin-once-removed label, preserving seven meioses. These pedigree categories are not biologically equivalent.

## Task 11 — Relatedness in inbred dog families

The corrected context model was a scikit-learn ExtraTrees classifier with 800 trees, minimum leaf size three, maximum feature fraction 0.75 and balanced sample weights. Its 185 features combined hap-IBD diploid/HBD summaries, chromosome-wise parental-haplotype concentration and 36 sibling-context summaries across 2/5/10 cM HBD cutoffs. Sibling context used supported cliques; average-linkage clustering supplied eight families. The frozen Campbell2016 sex-averaged canFam3.1 map covered chromosomes 1–38. Development crossed coalescent population size, MCAR missingness, ancestral mating and pedigree topology, correcting an earlier confounded schedule. The fixed comparison selected this candidate across 54 development cohorts: 89.46% exact accuracy, 80.81% macro recall and weighted Jaccard 0.8350, with temperature 0.75. Half-avuncular recall still failed strict gates; reserved stress/confirmation cohorts were not run. Final preprocessing cropped to map bounds and dropped incomplete sites without imputation; retained genotype/phase fields were preserved. All 68 samples entered eight families and all 265 inferred within-family pairs were exported. Finite templates, founder assumptions, true phase and MCAR missingness limit transfer.

## Task 12 — Relatedness with high background sharing

The fixed historical empirical-background histogram-gradient classifier used detected IBD summaries, 240 boosting iterations and the legacy background correction: a 25th-percentile reference, 5 cM minimum and corrected-feature weight one. Its original top-20 histogram truncation behavior was deliberately preserved because replacing it changed feature semantics; two replacement fits had failed to improve recorded development comparisons. For final preprocessing, 148,382 missing alleles were deterministically completed with the observed site-major allele, after preserving the original strict missingness failure. This unvalidated adaptation can bias IBD. All 43,956 pairs of 297 samples were retained. Following an unscored schema rejection, unsupported labels were merged into the unique accepted labels with the same meiosis count: 4,275 first-cousin-thrice-removed and 3,025 avuncular-twice-removed rows. This preserves meioses, not biological equivalence. The complete genotype-to-prediction pipeline has not passed independent scientific qualification.

## Task 13 — Low-coverage bottleneck

ANGSD genotype likelihoods (GL2) were integrated into site-allele-frequency likelihoods, including zero-depth uncertainty, rather than treating uncertain calls as exact genotypes. winsfs estimated the SFS; callable exposure and monomorphic mass were retained, and folding was global. moments fitted the two-epoch model with L-BFGS-B and an independent Powell check, scaling theta by mutation rate and callable length. Among existing fixed development methods, the continued winsfs run with 500 blocks, a 100-block window, tolerance 10⁻⁴ and a 50-epoch cap was selected before final input access. This selection used development error and numerical reliability, not final score feedback. The selected development fraction error was approximately +41%; a later full-batch refinement was worse and was not selected. This entry remains scientifically unqualified, particularly for ancestral-size/time and decline-fraction identifiability.

## Task 14 — Low-coverage split without migration

The same ANGSD GL2 likelihood approach produced a joint SFS, globally folded and fitted with a moments split-without-migration model. The fixed stabilized winsfs policy used 500 blocks and a full 500-block window, tolerance 10⁻⁶ and a 500-epoch cap. Two starts and continuation were retained to check spectrum and physical-target stability; the primary seed was fixed before final inference, with no bias correction or estimate selection on final outcomes. Both retained full-length development cohorts passed this policy's numerical and recovery checks, but the original six confirmation cohorts had failed an earlier policy, and no fresh six-cohort confirmation of the final policy was completed. Scientific qualification therefore remains false. An exact population-file header was removed in a separate audited staging step; BAM and reference identities were preserved.

## Data, reproducibility and limitations

Inputs comprised the official GHIST genotype, IBD and low-coverage files and their supplied references, callable masks, genetic maps and annotations. Development additionally used synthetic data from msprime/tskit, Ped-sim and stdpopsim/SLiM, with existing frozen human and dog maps and documented recombination/selection resources. Synthetic training and validation data are not independent empirical cohorts. The inherited Task 12 empirical-background profile was originally estimated from unlabeled challenge and reference cohorts; its complete historical provenance was not re-established in this review. This writeup does not assert that all development resources were challenge-provided or that no external data were used.

Run records retain source and model hashes, random seeds, input versions/checksums, optimizer failures, preprocessing changes and output audits. These records do not by themselves constitute a publicly reproducible release. A public code archive has not yet been verified; no public repository availability is claimed here. Some public testing-answer material was inadvertently displayed during the low-coverage specification review and was excluded from subsequent development. Consequently, the whole development history is not described as blinded. No final challenge truth was used for these methods.

The main lessons were that complete input exposure, correct coordinate and relationship-label conventions, numerical convergence and independent whole-cohort recovery need separate checks. Operational success did not prevent poor transfer. Explicit machine-readable task schemas, documented coordinate conventions and dashboard tables consistently named for their evaluation queues would reduce ambiguity in future rounds.

## Software references

- [moments](https://momentsld.github.io/moments/) and [dadi](https://dadi.readthedocs.io/en/latest/) for demographic SFS inference.
- [hap-IBD](https://github.com/browning-lab/hap-ibd) for phased-genotype IBD detection; [Ped-sim](https://github.com/williamslab/ped-sim) for pedigree simulation.
- [ANGSD](https://github.com/ANGSD/angsd) and [winsfs](https://github.com/malthesr/winsfs) for likelihood-based low-coverage spectrum estimation.
- [RAiSD](https://github.com/alachins/raisd) for sweep summaries; [msprime](https://tskit.dev/msprime/docs/stable/) and [stdpopsim](https://popsim-consortium.github.io/stdpopsim-docs/stable/) for population-genetic simulation.
- [PONDEROSA](https://doi.org/10.1093/genetics/iyaf094) motivated the haplotype-concentration feature; this was a custom implementation, not an execution of PONDEROSA.

## Author statement

Kevin Korfmann is the named participant. AI agents assisted with analysis, code development, execution, validation, independent review and drafting under his direction. This statement does not attribute agent-executed experiments or reviews to manual human work. A public code release is not currently claimed; the methods, recorded failures and preprocessing decisions are described above.
