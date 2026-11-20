"""Reject artifacts generated with the historical haploid population scale."""
from collections.abc import Mapping

import numpy as np


def require_diploid_scale(artifact: Mapping, *, context: str) -> None:
    value = np.asarray(artifact.get("coalescent_ploidy", 0))
    if value.shape != () or value.item() != 2:
        raise ValueError(
            f"{context} lacks verified diploid coalescent scaling. Regenerate "
            "simulations and retrain/recalibrate; legacy ploidy=1 artifacts "
            "cannot certify a diploid effective-size prior."
        )
