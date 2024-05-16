"""Useful utilities involving field of points."""

from collections.abc import Sequence

import numpy as np


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
