"""Detect pockets using a field of points."""

from typing import Any

from collections.abc import MutableMapping, Sequence

import MDAnalysis as mda
from loguru import logger

from ..fop.clean import remove_convex_fop, remove_fop_atom_clash
from ..fop.cluster import cluster_fop
from ..fop.compute import get_fop_inputs
from ..fop.gen import gen_fop
from ..utils.spatial import calculate_distance_two_points
from .utils import link_resnames


def get_fop_pocket(
    protein: Sequence[Sequence[float]],
    pocket_c_alphas: Sequence[Sequence[float]],
    center: Sequence[float],
    resolution: float,
    radius: float,
    clustering_model: type | None = None,
    clustering_kwargs: MutableMapping[str, Any] | None = None,
) -> Sequence[Sequence[float]]:
    r"""Calculate the field of points (FOP) for a pocket.

    This function generates a field of points (FOP) representing the shape
    of a pocket using the provided protein and alpha carbon coordinates.
    It applies several processing steps including removing steric clashes
    and clustering.

    Args:
        protein: XYZ coordinates of protein atoms.
        pocket_c_alphas: XYZ coordinates of all alpha carbons in the pocket.
        center: XYZ center of the pocket.
        resolution: Resolution for the FOP.
        radius: Radius of sphere for FOP.
        clustering_model: Scikit-learn clustering model class to use.
        clustering_kwargs: Keyword arguments for `clustering_model` initialization.

    Returns:
        Pocket shape as a field of points.

    Example:
        >>> protein_coords = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
        >>> alpha_carbons = [[1.1, 2.1, 3.1], [4.1, 5.1, 6.1]]
        >>> center = [2.0, 3.0, 4.0]
        >>> resolution = 0.5
        >>> radius = 5.0
        >>> fop = get_fop_pocket(protein_coords, alpha_carbons, center, resolution, radius)
        >>> print(f"Field of Points: {fop}")
    """
    # Generate the initial FOP and snapped center
    fop_pocket, fop_center = gen_fop(center, resolution, radius)

    # Trim protein atoms to those within the FOP sphere
    trimmed_coords_protein = [
        x for x in protein if calculate_distance_two_points(x, fop_center) < radius
    ]

    # Remove points in FOP that have steric clashes with protein
    fop_pocket = remove_fop_atom_clash(fop_pocket, trimmed_coords_protein)

    # Remove points outside the convex hull created by the alpha carbons in the pocket
    fop_pocket = remove_convex_fop(fop_pocket, pocket_c_alphas)

    # Cluster the FOP points if a clustering model is provided
    fop_pocket = cluster_fop(
        fop_pocket,
        clustering_model=clustering_model,
        clustering_kwargs=clustering_kwargs,
    )

    return fop_pocket


def get_fop_pocket_convenience(
    atoms_frame: mda.AtomGroup,
    subpex_config: "SubpexConfig",  # type: ignore # noqa: F821
    *args: Any,
    **kwargs: Any,
) -> Sequence[Sequence[float]]:
    r"""Convenience function for getting the FOP for a pocket.

    This function simplifies the process of calculating the field of points
    (FOP) for a pocket using a simulation frame and the SuPEx configuration.

    Args:
        atoms_frame: Trajectory frame to analyze.
        subpex_config: SuPEx configuration containing selection strings and other parameters.
        *args: Additional positional arguments.
        **kwargs: Additional keyword arguments.

    Returns:
        Pocket shape as a field of points.

    Example:
        >>> u = mda.Universe("topology.pdb", "trajectory.dcd")
        >>> frame = u.select_atoms("protein")
        >>> config = SubpexConfig(pocket=PocketConfig(selection_str="protein"))
        >>> fop = get_fop_pocket_convenience(frame, config)
        >>> print(f"Field of Points: {fop}")
    """
    inputs = get_fop_inputs(
        atoms_frame=atoms_frame,
        pocket_selection=subpex_config.pocket.selection_str,
        subpex_config=subpex_config,
    )
    return get_fop_pocket(**inputs)


def get_residues_in_pocket(
    u: mda.Universe, center: Sequence[float], radius: float, resnames: Sequence[str]
) -> mda.AtomGroup:
    r"""Get residues within the pocket.

    Args:
        u: MDAnalysis universe containing the protein.
        center: XYZ coordinates of the pocket center.
        radius: Radius of the pocket to be considered.
        resnames: Residue names for selection.

    Returns:
        mda.AtomGroup: AtomGroup of residues within the pocket.
    """
    selection_str = f"point {center[0]} {center[1]} {center[2]} {radius} and ({link_resnames(resnames, linker='or')})"
    return u.select_atoms(selection_str)


def get_pocket_selection(
    u: mda.Universe,
    center: Sequence[float],
    radius: float,
    water_dist: float = 6.7,
    water_resnames: Sequence[str] = ["WAT", "HOH", "H2O"],
    selection_append: str = "and (not name H*)",
) -> str:
    r"""Get the initial selection of the pocket based on the specified radius.

    This function identifies the initial selection of pocket residues based
    on a user-specified radius, and optionally considers distance from water molecules
    as an approximation of surface residues.

    Args:
        universe: MDA universe containing the protein.
        center: XYZ coordinates of the pocket center.
        radius: Radius of the pocket to be considered.
        water_dist: Constraint to use for proximity to calculate surface residues.
        water_resnames: Residue names for water molecules.
        selection_append: String to append to pocket selection. By default we exclude
            hydrogen atoms.

    Returns:
        Selection string of the pocket as used in MDAnalysis.

    Example:
        >>> u = mda.Universe("topology.pdb", "trajectory.dcd")
        >>> center = [10.0, 10.0, 10.0]
        >>> radius = 5.0
        >>> selection = get_pocket_selection(u, center, radius)
        >>> print(f"Pocket selection: {selection}")

    Note:
        If water molecules are found within `radius` of the pocket center, we will
        try to limit the pocket residues to those near the surface. We approximate
        Residues are included if they are within `water_dist` away from any water
        molecule.
    """
    residues_pocket = u.select_atoms(
        f"protein and point {center[0]} {center[1]} {center[2]} {radius}"
    )
    water_in_pocket = get_residues_in_pocket(u, center, radius, water_resnames)

    if len(water_in_pocket.residues) == 0:
        logger.warning("There are no water molecules to obtain distances")
        logger.warning("Some buried residues may end up in the pocket selection")
        pocket_residues = list(
            set(residue.resid for residue in residues_pocket.residues)
        )
        selection_pocket = " or ".join(f"resid {resid}" for resid in pocket_residues)
    else:
        list_close_water = []
        for water_residue in water_in_pocket.residues:
            water = u.select_atoms(f"resid {water_residue.resid}")
            for protein_residue in residues_pocket.residues:
                residue = u.select_atoms(f"resid {protein_residue.resid}")
                distance = calculate_distance_two_points(
                    water.center_of_geometry(), residue.center_of_geometry()
                )
                if (
                    distance < water_dist
                    and protein_residue.resid not in list_close_water
                ):
                    list_close_water.append(protein_residue.resid)
        selection_pocket = " or ".join(f"resid {resid}" for resid in list_close_water)

    selection_pocket += " " + selection_append
    return selection_pocket
