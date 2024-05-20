"""Compute properties of pockets."""

from typing import Any

from collections.abc import MutableMapping, Sequence

import MDAnalysis as mda
import numpy as np
import numpy.typing as npt
import scipy as sp

from ..fop.compute import get_fop_inputs, get_fop_volume
from ..pocket.detect import get_fop_pocket, get_fop_pocket_convenience
from ..utils.spatial import calculate_distance_two_points, get_centroid, get_rmsd


def get_pocket_volume_convenience(
    atoms_frame: mda.AtomGroup,
    subpex_config: "SubpexConfig",  # type: ignore # noqa: F821
    atoms_ref: mda.AtomGroup | None = None,
    *args: Any,
    **kwargs: Any
) -> float:
    """This function simplifies the calculation of the volume of a pocket
    from a given trajectory frame. It wraps around the lower-level
    [`get_fop_volume`][fop.compute.get_fop_volume] function,
    abstracting away some of the complexity and inputs.

    Args:
        atoms_frame: Trajectory frame to analyze.
        subpex_config: SuPEx configuration containing
            selection strings and other parameters.
        atoms_ref: Reference frame for comparison,
            defaults to None.
        *args: Additional positional arguments.
        **kwargs: Additional keyword arguments.

    Returns:
        The calculated volume of the pocket.

    Example:
        >>> u = mda.Universe("topology.pdb", "trajectory.dcd")
        >>> frame = u.select_atoms("protein")
        >>> config = SubpexConfig(pocket=PocketConfig(selection_str="protein"))
        >>> volume = get_pocket_volume_convenience(frame, config)
        >>> print(f"Pocket volume: {volume:.2f}")
    """
    inputs = get_fop_inputs(
        atoms_frame, subpex_config.pocket.selection_str, subpex_config, *args, **kwargs
    )
    inputs["fop"] = get_fop_pocket(**inputs)
    return get_fop_volume(**inputs)


def get_pocket_rmsd(
    atoms_ref: npt.NDArray[np.float64], atoms_frame: npt.NDArray[np.float64]
) -> float:
    """Calculates the root-mean-square deviation (RMSD) between two sets of atomic coordinates.

    This function computes the RMSD between a reference structure and a current
    frame, assuming each atom has an equal mass of one. RMSD is a measure of
    the average distance between atoms of superimposed proteins.

    Args:
        atoms_ref: Atomic coordinates of the reference structure.
        atoms_frame: Atomic coordinates of the current frame.

    Returns:
        The RMSD between the reference and current frame.

    Example:
        >>> atoms_ref = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        >>> atoms_frame = np.array([[1.1, 2.1, 3.1], [4.1, 5.1, 6.1]])
        >>> rmsd = get_pocket_rmsd(atoms_ref, atoms_frame)
        >>> print(f"RMSD: {rmsd:.2f}")
    """
    rmsd = get_rmsd(atoms_ref, atoms_frame)
    return float(rmsd)


def get_pocket_rmsd_convenience(
    atoms_frame: mda.AtomGroup,
    subpex_config: "SubpexConfig",  # type: ignore # noqa: F821
    atoms_ref: mda.AtomGroup | None = None,
    *args: Any,
    **kwargs: Any
) -> float:
    """A convenience wrapper around
    [`get_pocket_rmsd`][pocket.compute.get_pocket_rmsd] using a simulation frame.

    This function simplifies the calculation of the RMSD between a reference
    frame and a given trajectory frame using the subpex configuration.

    Args:
        atoms_frame: Trajectory frame to analyze.
        subpex_config: SuPEx configuration containing
            selection strings and other parameters.
        atoms_ref: Reference frame for comparison, must be provided if not None.
        *args: Additional positional arguments.
        **kwargs: Additional keyword arguments.

    Returns:
        The RMSD between the reference and current frame.

    Raises:
        ValueError: If `atoms_ref` is None.

    Example:
        >>> u = mda.Universe("topology.pdb", "trajectory.dcd")
        >>> frame = u.select_atoms("protein")
        >>> ref = u.select_atoms("protein").positions
        >>> config = SubpexConfig(pocket=PocketConfig(selection_str="protein"))
        >>> rmsd = get_pocket_rmsd_convenience(frame, config, ref)
        >>> print(f"Pocket RMSD: {rmsd:.2f}")
    """
    if atoms_ref is None:
        raise ValueError("atoms_ref cannot be None")
    rmsd = get_pocket_rmsd(atoms_ref.positions, atoms_frame.positions)
    return rmsd


def get_pocket_rog(fop: Sequence[Sequence[float]], *args: Any, **kwargs: Any) -> float:
    """Calculates the radius of gyration for a field of points (FOP).

    The radius of gyration is a measure of the compactness of a set of points,
    assuming a uniform mass distribution.

    Args:
        fop: The field of points defining the pocket shape.
        *args: Additional positional arguments.
        **kwargs: Additional keyword arguments.

    Returns:
        The radius of gyration of the pocket.

    Example:
        >>> fop = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
        >>> rog = get_pocket_rog(fop)
        >>> print(f"Radius of Gyration: {rog:.2f}")
    """
    mass = len(fop)

    if mass == 0:
        # Sometimes there are no points.
        return -9999

    centroid = get_centroid(fop)
    gyration_radius = 0.0
    for i in fop:
        gyration_radius += (calculate_distance_two_points(i, centroid)) ** 2

    return float(np.sqrt(gyration_radius / mass))


def get_pocket_rog_convenience(
    atoms_frame: mda.AtomGroup,
    subpex_config: "SubpexConfig",  # type: ignore # noqa: F821
    atoms_ref: mda.AtomGroup | None = None,
    *args: Any,
    **kwargs: Any
) -> float:
    """A convenience wrapper around
    [`get_pocket_rog`][pocket.compute.get_pocket_rog] using a simulation frame.

    This function simplifies the calculation of the radius of gyration for a
    pocket from a given trajectory frame using the subpex configuration.

    Args:
        atoms_frame: Trajectory frame to analyze.
        subpex_config: SuPEx configuration containing
            selection strings and other parameters.
        atoms_ref: Reference frame for comparison, must be provided if not None.
        *args: Additional positional arguments.
        **kwargs: Additional keyword arguments.

    Returns:
        The radius of gyration of the pocket.

    Examples:
        >>> u = mda.Universe("topology.pdb", "trajectory.dcd")
        >>> frame = u.select_atoms("protein")
        >>> config = SubpexConfig(pocket=PocketConfig(selection_str="protein"))
        >>> rog = get_pocket_rog_convenience(frame, config)
        >>> print(f"Radius of Gyration: {rog:.2f}")
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
    """Calculates the Jaccard distance between two fields of points (FoPs).

    The Jaccard distance measures the dissimilarity between two sets of points
    based on their intersection and union, considering a specific resolution
    for point proximity.

    Args:
        fop_ref: Reference field of points.
        fop_frame: Field of points from the current frame.
        resolution: Resolution for determining point proximity.

    Returns:
        The Jaccard distance between the reference and current field of points.

    Examples:
        >>> fop_ref = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
        >>> fop_frame = [[1.1, 2.1, 3.1], [4.1, 5.1, 6.1]]
        >>> jaccard = get_pocket_jaccard(fop_ref, fop_frame, 0.5)
        >>> print(f"Jaccard Distance: {jaccard:.2f}")
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
    subpex_config: "SubpexConfig",  # type: ignore # noqa: F821
    atoms_ref: mda.AtomGroup | None = None,
    *args: Any,
    **kwargs: Any
) -> float:
    """A convenience wrapper around
    [`get_pocket_jaccard`][pocket.compute.get_pocket_jaccard] using a simulation frame.

    This function simplifies the calculation of the Jaccard distance between
    a reference frame and a given trajectory frame using the subpex configuration.

    Args:
        atoms_frame: Trajectory frame to analyze.
        subpex_config: SuPEx configuration containing
            selection strings and other parameters.
        atoms_ref: Reference frame for comparison,
            must be provided if not None.
        *args: Additional positional arguments.
        **kwargs: Additional keyword arguments.

    Returns:
        The Jaccard distance between the reference and current field of points.

    Examples:
        >>> u = mda.Universe("topology.pdb", "trajectory.dcd")
        >>> frame = u.select_atoms("protein")
        >>> ref = u.select_atoms("protein").positions
        >>> config = SubpexConfig(pocket=PocketConfig(selection_str="protein"))
        >>> jaccard = get_pocket_jaccard_convenience(frame, config, ref)
        >>> print(f"Pocket Jaccard Distance: {jaccard:.2f}")
    """
    inputs: MutableMapping[str, Sequence[Sequence[float]] | float] = {}
    inputs["fop_ref"] = get_fop_pocket_convenience(
        atoms_frame=atoms_ref, subpex_config=subpex_config
    )
    inputs["fop_frame"] = get_fop_pocket_convenience(
        atoms_frame=atoms_frame, subpex_config=subpex_config
    )
    inputs["resolution"] = subpex_config.pocket.resolution
    return get_pocket_jaccard(**inputs)  # type: ignore
