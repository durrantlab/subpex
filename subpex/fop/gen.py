"""Generate field of points."""

from collections.abc import Sequence

import numpy as np

from ..utils.spatial import calculate_distance_two_points


def gen_fop(
    center: Sequence[float], resolution: float, radius: float, round_decimals: int = 3
) -> tuple[Sequence[Sequence[float]], Sequence[float]]:
    """Generate a spherical field of points.

    Snaps the provided center to the a center congruent to the resolution
    and radius given by the user. Then generates a spherical field of points
    (FOP) with radius = radius/2 with points spaced "resolution" apart.

    Args:
        center: XYZ coordinates for the center of the pocket.
        resolution: Distance between points which give resolution of the
            pocket volume calculator method.
        radius: Defines the radius of the sphere of the FOP.
        round_decimals: Number of decimals to round FOP and center to.

    Returns:
        field_of_points (list of lists): a spherical field of points
        center (list): center of the field of points
    """
    # This part of the code will get a new center, a snapped center.
    snapped_x = np.around(center[0] / resolution) * resolution
    snapped_y = np.around(center[1] / resolution) * resolution
    snapped_z = np.around(center[2] / resolution) * resolution
    fop_center = np.around(
        np.array([snapped_x, snapped_y, snapped_z], dtype=np.float64), round_decimals
    ).tolist()

    # Getting the field of points around the center.
    xs = np.round(
        np.arange(
            fop_center[0] - radius, fop_center[0] + radius + resolution, resolution
        ),
        round_decimals,
    )
    ys = np.round(
        np.arange(
            fop_center[1] - radius, fop_center[1] + radius + resolution, resolution
        ),
        round_decimals,
    )
    zs = np.round(
        np.arange(
            fop_center[2] - radius, fop_center[2] + radius + resolution, resolution
        ),
        round_decimals,
    )
    fop = np.vstack(np.meshgrid(xs, ys, zs)).reshape(3, -1).T.tolist()

    # making the field of points a sphere of radius = radius/2.
    fop = [x for x in fop if calculate_distance_two_points(x, fop_center) < radius]

    return fop, fop_center
