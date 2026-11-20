import hashlib
import json
from pathlib import Path
import numpy as np
from ghist_secondary_contact.vcf_fast import scan_secondary_contact_vcf
root = Path("/vast/projects/smathi/cohort/korfmann/ghist2026")
vcf = root / "data/raw/demography/secondary_contact/testing/secondary_contact.testing.vcf.gz"
archive = root / "artifacts/task2_secondary_contact/moments_generator_p40x30_single_refined_unfolded.sfs.npz"
spectrum, metadata = scan_secondary_contact_vcf(vcf)
with np.load(archive) as data:
    identical = bool(np.array_equal(spectrum, data["spectrum"]))
result = {"scope": "Archived testing input regression; not current final data", "spectrum_matches_archived_extractor": identical, "spectrum_sha256": hashlib.sha256(spectrum.tobytes()).hexdigest(), "population": metadata}
Path("artifacts/reliability/task2_scan.json").write_text(json.dumps(result, indent=2)+"\n")
print(json.dumps({k: v for k, v in result.items() if k != "population"}))
if not identical:
    raise RuntimeError("Spectrum differs from archived extractor")
