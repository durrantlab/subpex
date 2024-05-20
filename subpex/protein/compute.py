from typing import Any

import MDAnalysis as mda
from MDAnalysis.analysis import align

from ..utils.spatial import get_rmsd


def get_backbone_rmsd(
    atoms_ref: mda.AtomGroup,
    subpex_config: "SubpexConfig",  # noqa: F821
    atoms_frame: mda.AtomGroup,
    *args: Any,
    **kwargs: Any
) -> float:
    """
    Compute the Root Mean Square Deviation (RMSD) of the backbone atoms between
    two molecular dynamics frames after alignment.

    This function aligns the backbone atoms of the current frame to a reference frame
    and then calculates the RMSD between the aligned structures. It uses the MDAnalysis
    library for alignment and the `get_rmsd` function for RMSD calculation.

    Args:
        atoms_ref: An MDAnalysis AtomGroup object representing the reference structure.
        atoms_frame: An MDAnalysis AtomGroup object representing the
            current frame to be aligned and compared.

    Returns:
        float: The RMSD value between the backbone atoms of the reference and current frame.

    Example:
        >>> import MDAnalysis as mda
        >>> u = mda.Universe('topology.psf', 'trajectory.dcd')
        >>> ref = mda.Universe('reference.pdb')
        >>> rmsd = get_backbone_rmsd(ref, u)
        >>> print(f"Backbone RMSD: {rmsd:.2f}")
    """
    align.alignto(atoms_frame, atoms_ref, select="backbone")
    rmsd = get_rmsd(atoms_ref.positions, atoms_frame.positions)
    return float(rmsd)
