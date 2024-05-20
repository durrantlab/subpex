"""Generate field of points."""

from collections.abc import Sequence

import numpy as np
from loguru import logger


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
        field_of_points: a spherical field of points
        center: center of the field of points
    """
    logger.info("Generating spherical field of points (FOP)")
    logger.debug(f"FOP radius: {radius}")
    logger.debug(f"FOP resolution: {resolution}")
    # This part of the code will get a new center, a snapped center.
    snapped_x = np.around(center[0] / resolution) * resolution
    snapped_y = np.around(center[1] / resolution) * resolution
    snapped_z = np.around(center[2] / resolution) * resolution
    fop_center_array = np.around(
        np.array([snapped_x, snapped_y, snapped_z], dtype=np.float64), round_decimals
    )
    logger.debug(
        f"FOP center: {fop_center_array[0]} {fop_center_array[1]} {fop_center_array[2]}"
    )

    # Getting the field of points around the center.
    xs = np.round(
        np.arange(
            fop_center_array[0] - radius,
            fop_center_array[0] + radius + resolution,
            resolution,
        ),
        round_decimals,
    )
    ys = np.round(
        np.arange(
            fop_center_array[1] - radius,
            fop_center_array[1] + radius + resolution,
            resolution,
        ),
        round_decimals,
    )
    zs = np.round(
        np.arange(
            fop_center_array[2] - radius,
            fop_center_array[2] + radius + resolution,
            resolution,
        ),
        round_decimals,
    )
    fop_array = np.vstack(np.meshgrid(xs, ys, zs)).reshape(3, -1).T
    logger.debug(f"Initial FOP has {fop_array.shape[0]} points")
    dist_squared = np.sum((fop_array - fop_center_array) ** 2, axis=1)
    fop_array = fop_array[dist_squared < radius**2]
    logger.debug(f"{fop_array.shape[0]} points are < {radius} away from center")

    return fop_array.tolist(), fop_center_array.tolist()
