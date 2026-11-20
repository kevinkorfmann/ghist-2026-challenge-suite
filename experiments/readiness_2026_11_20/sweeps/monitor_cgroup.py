"""Read-only pilot memory sampling inside an already allocated Slurm node."""
import json
import sys
import time
from pathlib import Path

root = Path(f"/sys/fs/cgroup/memory/slurm/uid_1068187619/job_{sys.argv[1]}")
previous = None
for index in range(2400):
    try:
        values = {
            key: (root / f"memory.{key}").read_text().strip()
            for key in ("failcnt", "max_usage_in_bytes", "usage_in_bytes", "limit_in_bytes", "oom_control")
        }
    except FileNotFoundError:
        print(json.dumps({"unix_time": time.time(), "state": "cgroup_removed"}), flush=True)
        break
    if index % 10 == 0 or previous is None or values["oom_control"] != previous["oom_control"] or values["failcnt"] != previous["failcnt"]:
        print(json.dumps({"unix_time": time.time(), **values}), flush=True)
    previous = values
    time.sleep(1)
