"""Compute descriptors for pockets."""

import MDAnalysis as mda
import numpy as np
import numpy.typing as npt


def get_rmsd(R: npt.NDArray[np.float64], R_ref: npt.NDArray[np.float64]) -> float:
    """Compute atom RMSD.

    Args:
        R: Atomic coordinates.
        R_ref: Atomic coordinates of reference.
    """
    rmsd = mda.analysis.rms.rmsd(R_ref, R)
    return float(rmsd)
