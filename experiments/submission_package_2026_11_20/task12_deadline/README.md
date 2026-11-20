# Task 12 deadline candidate evidence

This archive records the fixed historical model used for the user-authorized best-effort final candidate. The input is current Synapse `syn75203676` v1. No training or biological parameter search was performed on the final data. Scientific qualification remains false.

The first strict attempt, Slurm 10859426, stopped on missing genotypes before detector/model inference. Its source, lock and failure evidence are retained under `original_strict/`. The v2 operational adaptation completes observed missing alleles with the observed per-site major allele: 148,382 of 147,944,610 alleles (approximately 0.1003%). Three multiallelic sites are skipped by the existing preparation route. Original observed alleles and phasing guards remain in place.

Job 10859506 completed in 2m28s with exit 0. Its candidate `v2/output/task12.final.csv` covers all 43,956 unordered pairs across 297 input samples. The independent schema/pair/hash check and original source manifests are retained. All pairs had detected segments, so the zero-segment fallback was unused. The historical correction retains its documented histogram truncation defect to preserve the selected model's feature semantics; its accuracy is not established for this final dataset.

`candidate_manifest.json` identifies the selected CSV by SHA256 and enumerates remaining scientific limitations. Model assets and final input VCFs are retained in the verified execution environment, with hashes and locations in the method/source records. The 39.9 MB probability artifact remains in safe local cache and is indexed in `large_artifact_index.json`; large prepared VCF and detector outputs remain on Betty and are bound by `v2/output/remote_output_manifest.json`.

These are executed-run records. The archived runner's relative paths refer to its original frozen execution directory; copying this evidence directory alone does not stage the external model/data assets. No competition upload or submission was performed by these runners. The central submission owner records receipts separately.
