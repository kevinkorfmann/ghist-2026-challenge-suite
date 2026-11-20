"""Version simulation features/calibration after scientific extraction fixes."""

SPECTRUM_CONTRACT = "joint_fixed_differences_retained_epochs_v2"


def require_spectrum_contract(record):
    if str(record.get("spectrum_contract", "")) != SPECTRUM_CONTRACT:
        raise ValueError(
            "Regenerate Task 2 simulation/calibration or retrain checkpoint: fixed-difference spectrum contract is missing or obsolete"
        )
    return SPECTRUM_CONTRACT
