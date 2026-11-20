"""Offline checks of the checked-in jobs; never execute or submit a job.

These checks cover script headers, not sbatch environment/CLI overrides or the
live partition defaults. Missing memory requests deliberately use those defaults.
"""

from pathlib import Path
import re
import shlex
import subprocess

import pytest


ROOT = Path(__file__).resolve().parents[1]
JOBS = sorted((ROOT / "slurm").glob("*.sbatch"))
MEMORY_MIB = {
    "genoa-std-mem": 5632,
    "genoa-lrg-mem": 15872,
    "b200-mig45": 8192,
    "b200-mig90": 8192,
    "dgx-b200": 8192,
}


def header_resources(text):
    """Read resource options in the single-component headers used here.

    The directive marker must start in column one. Slurm stops processing at
    the first executable line; quoted hashes are values, unquoted hashes are
    comments. This is intentionally not a shell parser or an sbatch validator.
    """
    options = {}
    names = {
        "--partition": "partition", "-p": "partition",
        "--mem-per-cpu": "mem-per-cpu", "--mem": "mem",
        "--mem-per-gpu": "mem-per-gpu",
    }
    # Consume the value of unrelated options too, so e.g. a --comment value
    # containing '--mem=...' is never interpreted as a resource request.
    valued = {
        "--job-name", "-J", "--cpus-per-task", "-c", "--time", "-t",
        "--output", "-o", "--error", "-e", "--array", "-a", "--gres",
        "--comment", "--chdir", "-D", "--account", "-A", "--nodes", "-N",
        "--ntasks", "-n", "--export", "--constraint", "-C", "--dependency",
        "-d", "--qos", "-q",
    } | names.keys()
    for line in text.splitlines():
        if line.startswith(("#SBATCH", "#SLURM")):
            prefix = "#SBATCH" if line.startswith("#SBATCH") else "#SLURM"
            tokens = iter(shlex.split(line[len(prefix):], comments=True))
            for token in tokens:
                option, separator, value = token.partition("=")
                if token.startswith("-p") and not token.startswith("--"):
                    option, value = "-p", token[2:]
                    separator = bool(value)
                if option in valued and not separator:
                    value = next(tokens, None)
                    assert value is not None, f"Missing value in {line!r}"
                if option in names:
                    options[names[option]] = value
        elif line.strip() and not line.lstrip().startswith("#"):
            break
    return options


def memory_mib(value):
    match = re.fullmatch(r"\s*(\+?\d+)(?:([KMGT])B?)?", value, re.IGNORECASE)
    assert match, f"Unsupported memory size {value!r}"
    size, unit = int(match[1]), (match[2] or "").upper()
    if unit == "K":
        return (size + 1023) // 1024
    return size * {"": 1, "M": 1, "G": 1024, "T": 1024**2}[unit]


def check_memory_policy(text):
    options = header_resources(text)
    requests = {"mem", "mem-per-cpu", "mem-per-gpu"} & options.keys()
    assert len(requests) <= 1, f"Mutually exclusive memory options: {sorted(requests)}"
    # Unspecified memory and per-node/per-GPU requests need live scheduler
    # configuration to assess; do not pretend that an offline check covers them.
    if "mem-per-cpu" in options:
        for partition in options.get("partition", "").split(","):
            if partition in MEMORY_MIB:
                assert memory_mib(options["mem-per-cpu"]) == MEMORY_MIB[partition], (
                    f"{partition} requires {MEMORY_MIB[partition]}M per CPU"
                )


def test_job_inventory_is_not_empty():
    assert JOBS, "No slurm/*.sbatch files found; refusing a vacuous audit"


@pytest.mark.parametrize("path", JOBS, ids=lambda path: path.name)
def test_memory_requests_follow_partition_policy(path):
    check_memory_policy(path.read_text())


@pytest.mark.parametrize("path", JOBS, ids=lambda path: path.name)
def test_jobs_have_valid_bash_syntax(path):
    # BASH_ENV can otherwise execute arbitrary startup code, even under bash -n.
    result = subprocess.run(
        ["bash", "--noprofile", "--norc", "-n", str(path)],
        env={"PATH": "/usr/bin:/bin"}, capture_output=True, text=True, check=False,
    )
    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize("partition", [
    "--partition=genoa-std-mem", "--partition genoa-std-mem",
    "-p genoa-std-mem", "-pgenoa-std-mem",
])
@pytest.mark.parametrize("memory", ["--mem-per-cpu=5500M", "--mem-per-cpu 5500M"])
def test_resource_formatting_cannot_bypass_policy(partition, memory):
    with pytest.raises(AssertionError, match="5632M"):
        check_memory_policy(f"#SBATCH {partition} {memory}\ntrue\n")


def test_header_boundary_quotes_and_inline_comments():
    text = '''#!/bin/bash
# Ordinary comment
#SBATCH --comment "--mem=10G # just text"
#SBATCH --partition 'b200-mig45' # partition comment
#SBATCH --mem-per-cpu="8G"# memory comment
    #SBATCH --mem=999G
set -euo pipefail
#SBATCH --mem-per-cpu=1M
python - <<'PY'
#SBATCH --mem=2G
PY
'''
    assert header_resources(text) == {"partition": "b200-mig45", "mem-per-cpu": "8G"}
    check_memory_policy(text)


@pytest.mark.parametrize("value", ["8192M", "8192", "8G", "8388608K", "8gb", "+8192M"])
def test_equivalent_memory_units(value):
    check_memory_policy(f"#SBATCH -p b200-mig45 --mem-per-cpu={value}\n")


@pytest.mark.parametrize("option", ["--mem=8G", "--mem-per-gpu 8G"])
def test_mutually_exclusive_memory_requests_are_rejected(option):
    with pytest.raises(AssertionError, match="Mutually exclusive"):
        check_memory_policy(f"#SBATCH --mem-per-cpu=8192M {option}\n")


def test_missing_memory_request_can_use_scheduler_defaults():
    check_memory_policy("#SBATCH --partition=genoa-std-mem\ntrue\n")


def test_shell_variables_in_directives_are_not_expanded():
    with pytest.raises(AssertionError, match="Unsupported memory size"):
        check_memory_policy("#SBATCH -p b200-mig45 --mem-per-cpu=$MEMORY\n")
