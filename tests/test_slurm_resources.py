"""Prevent reintroducing jobs rejected by Betty's enforced memory ratios."""

from pathlib import Path
import re


def test_memory_requests_follow_partition_policy():
    expected = {
        "genoa-std-mem": "5632M",
        "genoa-lrg-mem": "15872M",
        "b200-mig45": "8192M",
        "b200-mig90": "8192M",
        "dgx-b200": "8192M",
    }
    for path in Path("slurm").glob("*.sbatch"):
        text = path.read_text()
        partition = re.search(r"^#SBATCH --partition=(.+)$", text, re.M)
        memory = re.search(r"^#SBATCH --mem-per-cpu=(.+)$", text, re.M)
        if partition and partition[1] in expected and memory:
            assert memory[1] == expected[partition[1]], str(path)
