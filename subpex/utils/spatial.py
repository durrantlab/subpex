"""Useful utilities involving field of points."""

from collections.abc import Sequence

import numpy as np
import numpy.typing as npt


def get_centroid(fop_pocket: Sequence[Sequence[float]]) -> Sequence[float]:
    """Calculates the centroid or the center of mass for equal masses of a
    field of points.

    Args:
        fop_pocket:  the field of points defining the pocket shape.

    Returns:
        Coordinates of the pocket centroid.
    """
    x = float(np.mean([x[0] for x in fop_pocket]))
    y = float(np.mean([y[1] for y in fop_pocket]))
    z = float(np.mean([z[2] for z in fop_pocket]))
    return [x, y, z]


def calculate_distance_two_points(
    point1: Sequence[float], point2: Sequence[float]
) -> float:
    """Calculates the distance between two points in a 3D coordinate system.

    Args:
        point1: coordinates of first point. [x, y, z]
        point2: coordinates of second point. [x, y, z]

    Returns:
        Distance between `point1` and `point2`.
    """

    distance = 0.0
    for p1, p2 in zip(point1, point2):
        distance += (p1 - p2) ** 2

    return float(np.sqrt(distance))


def get_rmsd(R: npt.NDArray[np.float64], R_ref: npt.NDArray[np.float64]) -> float:
    """Compute the Root Mean Square Deviation (RMSD) between two sets of atomic coordinates.

    RMSD is a measure of the average distance between atoms of two superimposed
    sets of coordinates. This function calculates RMSD without any alignment or
    fitting of the structures, assuming that the input arrays are already aligned.

    Args:
        R (np.ndarray): A 2D array of shape (N, 3) representing the atomic coordinates.
                        Each row corresponds to the coordinates of one atom.
        R_ref (np.ndarray): A 2D array of shape (N, 3) representing the reference atomic coordinates.
                            Each row corresponds to the coordinates of one atom.

    Returns:
        float: The RMSD value between the two sets of atomic coordinates.

    Raises:
        ValueError: If the input arrays R and R_ref do not have the same shape.

    Example:
        >>> R = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        >>> R_ref = np.array([[1.1, 2.1, 3.1], [4.1, 5.1, 6.1]])
        >>> get_rmsd(R, R_ref)
        0.1
    """
    if R.shape != R_ref.shape:
        raise ValueError("R and R_ref must have the same shape")

    diff = R - R_ref
    diff_squared = np.square(diff)
    mean_diff_squared = np.mean(diff_squared)
    rmsd = np.sqrt(mean_diff_squared)

    return float(rmsd)
