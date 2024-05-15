"""Field of points properties"""

from collections.abc import Sequence


def get_fop_volume(
    fop: Sequence[Sequence[float]],
    resolution: float,
) -> float:
    """Compute field of point's volume.

    Args:
        fop: XYZ coordinates for each atom.
        resolution: resolution in Angstroms.
    """
    return float(len(fop) * (resolution**3))
