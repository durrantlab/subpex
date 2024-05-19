"""Compute properties of pockets."""

from typing import Any

from collections.abc import Sequence

import MDAnalysis as mda
import numpy as np
import numpy.typing as npt
import scipy as sp

from ..configs import SubpexConfig
from ..fop.compute import get_fop_inputs
from ..pocket.detect import get_fop_pocket
from ..utils.spatial import calculate_distance_two_points, get_centroid, get_rmsd


def get_pocket_rmsd(
    atoms_ref: npt.NDArray[np.float64], atoms_frame: npt.NDArray[np.float64]
) -> float:
    """Takes the xyz coordinates of a field of points and calculates the radius
    of gyration. It assumes a mass of one for all the points.

    Args:
        atoms_ref npt.NDArray[np.float64] : Atomic coordinates of reference structure.
        atoms_frame npt.NDArray[np.float64] : Atomic coordinates of current frame.

    Returns:
        RMSD of atoms captured by the pocket.
    """
    rmsd = get_rmsd(atoms_ref, atoms_frame)
    return float(rmsd)


def get_pocket_rmsd_convenience(
    atoms_frame: mda.AtomGroup,
    subpex_config: SubpexConfig,
    atoms_ref: mda.AtomGroup | None = None,
    *args: Any,
    **kwargs: Any
) -> float:
    """A convenience wrapper around `get_pocket_rmsd` using a simulation frame.

    Args:
        atoms_frame: Trajectory frame to analyze.
        subpex_config: SuPEx configuration.
    """
    if atoms_ref is None:
        raise ValueError("atoms_ref cannot be None")
    rmsd = get_pocket_rmsd(atoms_ref, atoms_frame)
    return rmsd


def get_pocket_rog(
    fop_pocket: Sequence[Sequence[float]], *args: Any, **kwargs: Any
) -> float:
    """Takes the xyz coordinates of a field of points and calculates the radius
    of gyration. It assumes a mass of one for all the points.

    Args:
        fop_pocket Sequence[Sequence[float]] : The field of points defining the pocket shape.

    Returns:
        radius of gyration of the pocket.
    """
    mass = len(fop_pocket)

    if mass == 0:
        # Sometimes there are no points.
        return -9999

    centroid = get_centroid(fop_pocket)
    gyration_radius = 0.0
    for i in fop_pocket:
        gyration_radius += (calculate_distance_two_points(i, centroid)) ** 2

    return float(np.sqrt(gyration_radius / mass))


def get_pocket_rog_convenience(
    atoms_frame: mda.AtomGroup,
    subpex_config: SubpexConfig,
    atoms_ref: mda.AtomGroup | None = None,
    *args: Any,
    **kwargs: Any
) -> float:
    """A convenience wrapper around `get_pocket_rog` using a simulation frame.

    Args:
        atoms_frame: Trajectory frame to analyze.
        subpex_config: SuPEx configuration.
    """
    inputs = get_fop_inputs(
        atoms_frame, subpex_config.pocket.selection_str, subpex_config, *args, **kwargs
    )
    inputs["fop"] = get_fop_pocket(**inputs)
    return get_pocket_rog(**inputs)


def get_pocket_jaccard(
    fop_ref: Sequence[Sequence[float]],
    fop_frame: Sequence[Sequence[float]],
    resolution: float,
) -> float:
    """Calculates the Jaccard distance between the points in `fop_ref` and
    the `fop_segment`. Uses the distance between points to calculate the intersection.

    Args:
        fop_ref: Reference FOP.
        fop_frame: Segment FOP.
        resolution: Field of points resolution.

    Returns:
        Jaccard distance.
    """
    # sometimes no points are present in the FOP.
    if len(fop_frame) == 0:
        return 1.0

    # Obtaining the trees for both field of points
    tree_ref = sp.spatial.cKDTree(fop_ref)
    tree_segment = sp.spatial.cKDTree(fop_frame)

    # Obtain the points that are at less than resolution/2.5 (aka have the same
    # coordinates)
    clash_indices = tree_ref.query_ball_tree(tree_segment, resolution / 2.5, p=2, eps=0)

    # Count the points that intersect and convert to float
    intersection = float(len([x for x in clash_indices if x]))

    # Obtain the union of both FOP
    union = float(len(fop_ref) + len(fop_frame) - intersection)

    # Calculate Jaccard distance
    jaccard = 1.0 - intersection / union

    return jaccard


def get_pocket_jaccard_convenience(
    atoms_frame: mda.AtomGroup,
    subpex_config: SubpexConfig,
    atoms_ref: mda.AtomGroup | None = None,
    fop_ref: Sequence[Sequence[float]] | None = None,
    *args: Any,
    **kwargs: Any
) -> float:
    """A convenience wrapper around `get_pocket_rog` using a simulation frame.

    Args:
        atoms_frame: Trajectory frame to analyze.
        subpex_config: SuPEx configuration.
    """
    if fop_ref is None:
        raise ValueError("fop_ref cannot be None")
    inputs = get_fop_inputs(
        atoms_frame, subpex_config.pocket.selection_str, subpex_config, *args, **kwargs
    )
    inputs["fop_frame"] = get_fop_pocket(**inputs)
    return get_pocket_jaccard(fop_ref=fop_ref, **inputs)
