from pathlib import Path


def test_clean_spectrum_holdout_uses_reserved_replicates():
    text = Path("slurm/task10_clean_spectrum_holdout.sbatch").read_text()
    for replicate in range(16, 20):
        assert f"replicate_{replicate}.seg" in text
