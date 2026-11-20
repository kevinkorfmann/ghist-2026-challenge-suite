"""Trace numerical calls for one previously censored independent configuration."""

import faulthandler
import json
from pathlib import Path
import sys
import time

import scripts.validate_demography_independently as benchmark
import ghist_secondary_contact.spectrum as secondary
import ghist_split_migration_bgs.spectrum as split

faulthandler.enable()
faulthandler.dump_traceback_later(30, repeat=True)
trace = Path(sys.argv[1]).open("w", buffering=1)
for module, names in (
    (
        secondary,
        ("secondary_contact_spectrum", "generator_matched_secondary_contact_spectrum"),
    ),
    (split, ("isolation_with_migration_spectrum",)),
):
    for name in names:
        original = getattr(module, name)

        def traced(parameters, sample_sizes, _original=original, _name=name, **kwargs):
            started = time.monotonic()
            trace.write(
                json.dumps(
                    {
                        "event": "start",
                        "function": _name,
                        "parameters": list(parameters),
                        "sample_sizes": list(map(int, sample_sizes)),
                        "kwargs": kwargs,
                        "time": started,
                    }
                )
                + "\n"
            )
            value = _original(parameters, sample_sizes, **kwargs)
            trace.write(
                json.dumps(
                    {
                        "event": "finish",
                        "function": _name,
                        "elapsed": time.monotonic() - started,
                    }
                )
                + "\n"
            )
            return value

        setattr(module, name, traced)
config = benchmark.configurations()[int(sys.argv[2])]
print(benchmark.run_one((config, 0, sys.argv[3], 2000000, 12)), flush=True)
