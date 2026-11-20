# Tasks 5–8: submitted and scored — 20 November 2026

All four current-final-data inference runs completed, with independent input/source/output integrity and BED geometry checks. **All four entries were submitted by the central uploader and are SCORED.** Central uploader: **Check submission status**, chat `01a1f9e8-d561-7de1-98f1-0255e41ead65`. Candidate handoff completed at 2026-11-20T21:30:25.432591+00:00. The deadline is 23:59 UTC. Server-score receipts verified at 2026-11-20T21:52:19.529894+00:00.

Files: `experiments/final_submission_2026_11_20/sweeps_takeover/ready/`; exact paths, checksums, inputs and audit paths are in its `manifest.json` and `SHA256SUMS`.

| Task | Candidate | Final queue | Chr9 zero-based half-open intervals | Compute |
|---|---|---:|---|---|
| 5 | `task5_final.bed` | 9619656 | 104,670,000–104,790,000 | kkor@128.91.222.212, 19.77 seconds |
| 6 | `task6_final.bed` | 9619680 | 115,437,265–115,587,265 | betty, 32 seconds |
| 7 | `task7_final.bed` | 9619658 | 2,803,500–2,963,500 | betty, 153 seconds |
| 8 | `task8_final.bed` | 9619662 | 37,720,879–37,840,879; 83,003,141–83,123,141 | betty, 50 seconds |

## Integrity and execution

All four VCFs are current Synapse version2 and their MD5s match the signed-in file sidebar; current entity/version/filehandle metadata and local SHA256 bindings are retained. The shared accessibility mask is currentv1. Frozen inference source matches the pre-input method lock, with no new biological simulations or tuning on final outcomes. Independent checks verified every bound artifact on its execution host, the unchanged source and input hashes, chr9/contig/callable bounds, interval counts/non-overlap and geometry reconstructed from the retained raw predictions. The production validator passed all four candidates.

Task5 ran on the user-provided ARM64 workstation;12 frozen runner tests passed there. Its first detached launcher failed before starting inference due to a PID-file race; a corrected launcher ran once and completed exit0. Peak memory424244KiB. RAiSD is not invoked by this Task5 method; its common-wrapper binary fingerprint does not imply an ARM RAiSD execution.

Task6 job10859516 completed exit0. Task7 job10859524 failed because the launcher omitted the existing GSL library directory. Retry10859617 added that runtime lookup only and completed exit0; source, binary, input and method settings were unchanged. Task8 job10859658 used the corrected completed Task7 reference and completed exit0. Runtime records, original failure, shell/SBatch commands, source freeze and library hashes are preserved. Large raw VCFs/RAiSD reports remain on their original compute hosts; their full hashes were verified and retained.

## Scientific limitations to retain in submission methods

These are best-effort entries, **not scientifically qualified releases**. Task5's three-estimator consensus test fails on its final input (span165783.5bp versus50000bp criterion). Task6 uses the previously selected50-SNP H12 development method without independent confirmation. Task7 lacks matched biological recovery and has disagreement from its500kb composite estimator; the frozen ten-estimator median is preserved. Task8 uses an uncalibrated threshold transferred from Task7 and produced two intervals without a fallback. The official exon reference is a two-column Ensembl annotation; its standard1-based inclusive interpretation was converted to BED and recorded in `exon_normalization.json`. The challenge wiki does not explicitly declare that exon convention, so one-base boundary uncertainty remains. No final truth or published answers were used.

## Submission receipts

The central uploader submitted all four unchanged candidates. This chat verified that the candidate bytes match the central upload staging files and archived the server-score receipt records under `experiments/final_submission_2026_11_20/sweeps_takeover/submission_receipts/`. The ready manifest now binds each receipt and score to the original candidate hash. No duplicate submissions were made.

| Task | Submission ID | Server status | Recall | F1 |
|---|---:|---|---:|---:|
| 5 | 9785369 | SCORED | 0 | 0 |
| 6 | 9785370 | SCORED | 0 | 0 |
| 7 | 9785371 | SCORED | 1 | 0.555556 |
| 8 | 9785372 | SCORED | 0.060606 | 0.106952 |

Acceptance is distinct from predictive quality. Tasks 5 and 6 scored zero F1; all original scientific limitations above remain. Do not alter predictions using final results or submit duplicates.
