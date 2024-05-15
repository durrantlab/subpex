"""Algorithms to clean a field of points."""

from collections.abc import Sequence

import numpy as np
import scipy as sp


def remove_fop_atom_clash(
    fop: Sequence[Sequence[float]],
    atoms_points: Sequence[Sequence[float]],
    vdw_radius: float = 2.6,
) -> Sequence[Sequence[float]]:
    """Removes points that clash with atoms within the field of points.

    Args:
        fop: Field of points XYZ coordinates.
        atoms_points (list of lists): coordinates of atoms within the field of points.
        vdw_radius: Van der Waals radius for atoms to detect clashes. Defaults to
            the van der Waals radius of water and that of hydrogen atom.

    Returns:
        A field of points without atom clashes.
    """
    # Generating cKDTrees to obtain distances
    atoms_tree = sp.spatial.cKDTree(atoms_points)
    fop_tree = sp.spatial.cKDTree(fop)

    # get indices of fop_tree that are close to atoms_tree
    clash_indices = atoms_tree.query_ball_tree(fop_tree, vdw_radius, p=2, eps=0)

    # Creating the FOP that does not contain points that clash with the atoms
    # provided to the function.
    fop_trimmed = []
    for i, point in enumerate(fop):
        if i not in clash_indices:
            fop_trimmed.append(point)

    return fop_trimmed


def point_in_hull(
    point: Sequence[float], hull: sp.spatial.ConvexHull, tolerance: float = 1e-12
) -> bool:
    """Determines if a point is contained in a hull. A point is in the hull if
    and only if for every equation (describing the facets), the dot product
    between the point and the normal vector (eq[:-1]) plus the offset (eq[-1])
    is less than or equal to zero. You may want to compare to a small, positive
    constant tolerance = 1e-12 rather than to zero because of issues of
    numerical precision (otherwise, you may find that a vertex of the convex
    hull is not in the convex hull).

    Args:
        point: coordinates of the point to be considered.
        hull: Convex hulls in N dimensions.
        tolerance: tolerance for the calculation.

    Returns:
        Returns True if point in hull, False otherwise
    """
    return all((np.dot(eq[:-1], point) + eq[-1] <= tolerance) for eq in hull.equations)


def remove_convex_fop(
    fop_trimmed: Sequence[Sequence[float]], trimmed_alpha: Sequence[Sequence[float]]
) -> Sequence[Sequence[float]]:
    """Uses a convex hull to trim points from the field of points that are
    outside of the protein pocket.

    Args:
        fop: Field of points XYZ coordinates with no atom clashes.
        trimmed_alpha: alpha atoms to define the convex hull.

    Returns:
        points_in_hull: FOP trimmed of points outside the convex hull defined by the
            protein.
    """
    trimmed_alpha_convex = sp.spatial.ConvexHull(trimmed_alpha)
    points_in_hull = []
    for point in fop_trimmed:
        if point_in_hull(point, trimmed_alpha_convex):
            points_in_hull.append(point)
    return points_in_hull
