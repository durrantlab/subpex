"""Detect pockets using a field of points."""

from typing import Any

from collections.abc import MutableMapping, Sequence

import MDAnalysis as mda
from loguru import logger

from ..fop.clean import remove_convex_fop, remove_fop_atom_clash
from ..fop.cluster import cluster_fop
from ..fop.gen import gen_fop
from ..utils.spatial import calculate_distance_two_points


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


def get_pocket_selection(
    universe: mda.Universe,
    center: Sequence[float],
    radius: float,
    distance_constraint: float = 6.7,
) -> str:
    """Takes the protein and center of the pocket and gets the initial
    selection of the pockets based on the user-specified radius. If water
    molecules are present, it uses them as a first approximation of the surface
    residues.

    Args:
        universe: MDA universe with protein.
        center: XYZ coordinates of the pocket center.
        radius: radius of the pocket to be considered.
        distance_constraint: Constraint to use for proximity to calculate surface
            residues.

    Returns:
        Selection string of the pocket as used in MDAnalysis.
    """
    water_pocket = universe.select_atoms(
        f"point {center[0]} {center[1]} {center[2]} {radius} and resname WAT"
    )
    residues_pocket = universe.select_atoms(
        f"protein and point {center[0]} {center[1]} {center[2]} {radius}"
    )

    if len(water_pocket.residues) == 0:
        logger.warning(
            "There are no water molecules to obtain distances, some buried residues may end up in the pocket selection"
        )
        # obtain all residues in the pocket and create a selection pocket string.
        pocket_residues = []
        for residue in residues_pocket.residues:
            if residue.resid not in pocket_residues:
                pocket_residues.append(residue.resid)
        selection_pocket = f"resid {pocket_residues[0]} "
        for i in pocket_residues[1:]:
            selection_pocket += f" or resid {str(i)} "

    else:
        distances = []
        list_close_water = []
        # we want to make sure that the residues we are selecting are close to the
        # surface, we will approximate the surface accessible residues
        # using water molecules and a distance cutoff for that purpose.
        for i in water_pocket.residues:
            water = universe.select_atoms(f"resid {i.resid}")
            for j in residues_pocket.residues:
                residue = universe.select_atoms(f"resid {j.resid}")
                distance = calculate_distance_two_points(
                    water.center_of_geometry(), residue.center_of_geometry()
                )
                distances.append(distance)
                if distance < distance_constraint and j.resid not in list_close_water:
                    list_close_water.append(j.resid)

        # creating selection pocket string compatible with
        selection_pocket = f"resid {list_close_water[0]}"
        for i in list_close_water[1:]:
            selection_pocket += f" or resid {str(i)}"

            # save the pocket selection string so it can be used in other places
    selection_pocket += " and (not name H*)"

    return selection_pocket
