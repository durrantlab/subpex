"""Compute properties of pockets."""

from collections.abc import Sequence

import numpy as np

from ..utils.spatial import calculate_distance_two_points, get_centroid


def calculate_pocket_gyration(fop_pocket: Sequence[Sequence[float]]) -> float:
    """Takes the xyz coordinates of a field of points and calculates the radius
    of gyration. It assumes a mass of one for all the points.

    Args:
        fop_pocket: the field of points defining the pocket shape.

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
