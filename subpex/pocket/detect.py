"""Detect pockets using a field of points."""

from typing import Any

from collections.abc import MutableMapping, Sequence

from ..fop.clean import remove_convex_fop, remove_fop_atom_clash
from ..fop.cluster import cluster_fop
from ..fop.gen import gen_fop
from ..fop.utils import calculate_distance_two_points


def get_fop_pocket(
    protein: Sequence[Sequence[float]],
    alphas: Sequence[Sequence[float]],
    center: Sequence[float],
    resolution: float,
    radius: float,
    clustering_model: type | None = None,
    clustering_kwargs: MutableMapping[str, Any] | None = None,
) -> Sequence[Sequence[float]]:
    """Takes the coordinates of the protein and the coordinates of alpha
    carbons and calculates the field of points (FOP).

    Args:
        protein: XYZ coordinates of protein atoms.
        alphas: XYZ coordinates of all alpha carbons.
        center: XYZ center of the pocket.
        resolution: resolution for the FOP.
        radius: radius of sphere for FOP.
        clustering_model: Scikit-learn clustering model class to use.
        clustering_kwargs: Keyword arguments for `clustering_model` initialization.

    Returns:
        Pocket shape as a field of points.
    """
    # get starting FOP and snapped center
    fop_pocket, fop_center = gen_fop(center, resolution, radius)
    # trim protein atoms to those within the FOP sphere.
    trimmed_coords_protein = [
        x for x in protein if calculate_distance_two_points(x, fop_center) < radius
    ]
    # remove points in FOP that have steric clashes with protein.
    fop_pocket = remove_fop_atom_clash(fop_pocket, trimmed_coords_protein)
    # remove points outside the convex hull created by the alpha carbons in the pocket.
    fop_pocket = remove_convex_fop(fop_pocket, alphas)
    fop_pocket = cluster_fop(
        fop_pocket,
        clustering_model=clustering_model,
        clustering_kwargs=clustering_kwargs,
    )

    return fop_pocket
