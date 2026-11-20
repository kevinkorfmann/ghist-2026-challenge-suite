"""Record native adaptive non-progress and interrupt only the reproducer."""

import json
from pathlib import Path
import sys
import scripts.validate_demography_independently as benchmark
import ghist_secondary_contact.spectrum as secondary
import ghist_split_migration_bgs.spectrum as split

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
            if kwargs.get("adapt_dt"):
                repeats = {}
                trace.write(
                    json.dumps(
                        {
                            "event": "adaptive_start",
                            "function": _name,
                            "parameters": list(parameters),
                            "sample_sizes": list(map(int, sample_sizes)),
                            "kwargs": kwargs,
                        }
                    )
                    + "\n"
                )

                def watch(frame, event, arg):
                    if frame.f_code.co_name != "integrate_nD":
                        return None
                    if (
                        event == "line"
                        and frame.f_lineno == 1062
                        and frame.f_locals.get("neg")
                    ):
                        state = frame.f_locals
                        key = (state["t"], state["dt"])
                        repeats[key] = repeats.get(key, 0) + 1
                        if repeats[key] >= 16:
                            trace.write(
                                json.dumps(
                                    {
                                        "event": "native_nonprogress_confirmed",
                                        "time": float(state["t"]),
                                        "dt": float(state["dt"]),
                                        "target": float(state["Tmax"]),
                                        "repetitions": repeats[key],
                                    }
                                )
                                + "\n"
                            )
                            raise KeyboardInterrupt(
                                "Native adaptive integrator repeated identical t and dt 16 times"
                            )
                    return watch

                sys.settrace(watch)
            try:
                return _original(parameters, sample_sizes, **kwargs)
            finally:
                sys.settrace(None)

        setattr(module, name, traced)
config = benchmark.configurations()[int(sys.argv[2])]
print(benchmark.run_one((config, 0, sys.argv[3], 2000000, 12)), flush=True)
