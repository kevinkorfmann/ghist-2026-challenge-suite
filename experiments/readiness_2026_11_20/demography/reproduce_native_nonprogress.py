"""Safely reproduce moments' native adaptive loop independently of our wrappers."""

import hashlib
import inspect
import json
from pathlib import Path
import sys
import time

import moments
import numpy as np


class NonProgress(RuntimeError):
    pass


def reproduce(parameters, *, secondary_contact):
    nu1, nu2 = parameters[:2]
    fs = moments.Spectrum(moments.LinearSystem_1D.steady_state_1D(36))
    fs = moments.Manips.split_1D_to_2D(fs, 20, 16)
    if secondary_contact:
        isolation, duration, migration = parameters[2:]
        fs.integrate([nu1, nu2], isolation, dt_fac=0.02, adapt_dt=False)
    else:
        duration, migration = parameters[2:]
    function = moments.Integration.integrate_nD
    lines, first_line = inspect.getsourcelines(function)
    loop_line = next(
        first_line + index
        for index, line in enumerate(lines)
        if line.strip() == "while t < Tmax:"
    )
    repeated = {}
    record = {
        "parameters": parameters,
        "sample_sizes": [20, 16],
        "dt_fac": 0.02,
        "secondary_contact": secondary_contact,
    }

    def watch(frame, event, arg):
        if frame.f_code is not function.__code__:
            return None
        if (
            event == "line"
            and frame.f_lineno == loop_line
            and frame.f_locals.get("neg")
        ):
            state = frame.f_locals
            key = (state["t"], state["dt"])
            repeated[key] = repeated.get(key, 0) + 1
            if repeated[key] == 16:
                record.update(
                    status="native_nonprogress_confirmed",
                    time=float(state["t"]),
                    dt=float(state["dt"]),
                    target=float(state["Tmax"]),
                    identical_repeats=16,
                )
                raise NonProgress("Repeated the same rejected step 16 times")
        return watch

    started = time.perf_counter()
    sys.settrace(watch)
    try:
        fs.integrate(
            [nu1, nu2],
            duration,
            m=np.array([[0, migration], [migration, 0]]),
            dt_fac=0.02,
            adapt_dt=True,
        )
        record["status"] = "integration_completed"
    except NonProgress:
        pass
    finally:
        sys.settrace(None)
    record["elapsed_seconds"] = time.perf_counter() - started
    return record


records = [
    reproduce(
        [
            2.3195170721030354,
            63.313519790087035,
            0.00255781614605999,
            0.7031593982045962,
            37.616984248152015,
        ],
        secondary_contact=True,
    ),
    reproduce(
        [
            37.002346797132866,
            6.218138603544293,
            0.059115084943017114,
            22.980411155782058,
        ],
        secondary_contact=False,
    ),
]
source = Path(inspect.getsourcefile(moments.Integration.integrate_nD))
result = {
    "moments_version": moments.__version__,
    "native_integration_source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
    "records": records,
    "method": "A Python trace observes repeated t/dt at the native while-loop and interrupts after 16 identical rejected steps; no third-party source is changed.",
}
Path(__file__).with_name("native_nonprogress.json").write_text(
    json.dumps(result, indent=2) + "\n"
)
print(json.dumps(result, indent=2))
