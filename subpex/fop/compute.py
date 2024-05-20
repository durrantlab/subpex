"""Field of points properties"""

from typing import Any

from collections.abc import MutableMapping, Sequence

import MDAnalysis as mda


def get_fop_inputs(
    atoms_frame: mda.AtomGroup,
    pocket_selection: str | None,
    subpex_config: "SubpexConfig",  # type: ignore # noqa: F821
    *args: Any,
    **kwargs: Any
) -> MutableMapping[str, Any]:
    """
    Generates a dictionary of inputs for field of points (FoP) calculations.

    This function selects atoms from the provided MDAnalysis AtomGroup according to
    the given pocket selection string and the Subpex configuration. It returns a
    dictionary containing the positions of protein atoms, pocket C-alpha atoms,
    pocket radius, resolution, and pocket center.

    Args:
        atoms_frame: An MDAnalysis AtomGroup object containing the atoms of the
            molecular structure.
        pocket_selection: A string specifying the selection criteria for pocket atoms.
            This string should follow the MDAnalysis selection syntax.
        subpex_config: The Subpex configuration object containing parameters for
            pocket radius, resolution, and center.
        *args: Additional positional arguments.
        **kwargs: Additional keyword arguments.

    Returns:
        A dictionary with the following keys:

            - `protein`: Positions of protein atoms.
            - `pocket_c_alphas`: Positions of pocket C-alpha atoms.
            - `radius`: Radius of the pocket.
            - `resolution`: Resolution of the pocket.
            - `center`: Center of the pocket.

    Raises:
        ValueError: If the pocket selection string is `None` or if the pocket center
            is `None`.

    Example:
        ```python
        atoms_frame = u.select_atoms("protein")
        pocket_selection = "resid 50-60"
        subpex_config = SubpexConfig(...)
        inputs = get_fop_inputs(atoms_frame, pocket_selection, subpex_config)
        ```
    """
    if pocket_selection is None:
        raise ValueError("Pocket selection string cannot be `None`")
    inputs = {}
    inputs["protein"] = atoms_frame.select_atoms("protein").positions
    inputs["pocket_c_alphas"] = atoms_frame.select_atoms(
        pocket_selection + " and name CA*"
    ).positions
    inputs["radius"] = subpex_config.pocket.radius
    inputs["resolution"] = subpex_config.pocket.resolution
    if subpex_config.pocket.center is None:
        raise ValueError("`pocket_center` cannot be None")
    inputs["center"] = subpex_config.pocket.center
    return inputs


def get_fop_volume(
    fop: Sequence[Sequence[float]], resolution: float, *args: Any, **kwargs: Any
) -> float:
    """Computes the volume of a field of points (FoP) based on the given resolution.

    This function calculates the volume occupied by a field of points by assuming that
    each point occupies a cubic space defined by the resolution.

    Args:
        fop: A sequence of sequences containing XYZ coordinates for each point.
        resolution: The resolution of the field of points in Angstroms.
        *args: Additional positional arguments.
        **kwargs: Additional keyword arguments.

    Returns:
        The volume of the field of points as a float.

    Example:
        ```python
        fop = [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [2.0, 2.0, 2.0]]
        resolution = 1.0
        volume = get_fop_volume(fop, resolution)
        ```
    """
    return float(len(fop) * (resolution**3))
