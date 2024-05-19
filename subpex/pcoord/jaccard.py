from collections.abc import Sequence

import scipy as sp

from ..configs import SubpexConfig


def get_jaccard_distance(
    fop_ref: Sequence[Sequence[float]],
    fop_segment: Sequence[Sequence[float]],
    subpex_config: SubpexConfig,
) -> float:
    """Calculates the Jaccard distance between the points in `fop_ref` and
    the `fop_segment`. Uses the distance between points to calculate the intersection.

    Args:
        fop_ref: Reference FOP.
        fop_segment: Segment FOP.
        subpex_config: The subpex context manager.

    Returns:
        Jaccard distance.
    """
    # sometimes no points are present in the FOP.
    if len(fop_segment) == 0:
        return 1.0

    # Obtaining the trees for both field of points
    tree_ref = sp.spatial.cKDTree(fop_ref)
    tree_segment = sp.spatial.cKDTree(fop_segment)

    # Obtain the points that are at less than resolution/2.5 (aka have the same
    # coordinates)
    clash_indices = tree_ref.query_ball_tree(
        tree_segment, subpex_config.pocket_resolution / 2.5, p=2, eps=0
    )

    # Count the points that intersect and convert to float
    intersection = float(len([x for x in clash_indices if x]))

    # Obtain the union of both FOP
    union = float(len(fop_ref) + len(fop_segment) - intersection)

    # Calculate Jaccard distance
    jaccard = 1.0 - intersection / union

    return jaccard
