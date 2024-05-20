"""Field of points properties"""

from typing import Any

from collections.abc import MutableMapping, Sequence

import MDAnalysis as mda


def get_fop_inputs(
    atoms_frame: mda.AtomGroup,
    pocket_selection: str | None,
    subpex_config: "SubpexConfig",  # noqa: F821
    *args: Any,
    **kwargs: Any
) -> MutableMapping[str, Any]:
    if pocket_selection is None:
        raise ValueError("Pocket selection string cannot be `None`")
    data = {}
    data["protein"] = atoms_frame.select_atoms("protein").positions
    data["pocket_c_alphas"] = atoms_frame.select_atoms(
        pocket_selection + " and name CA*"
    ).positions
    data["radius"] = subpex_config.pocket.radius
    data["resolution"] = subpex_config.pocket.resolution
    if subpex_config.pocket.center is None:
        raise ValueError("`pocket_center` cannot be None")
    data["center"] = subpex_config.pocket.center
    return data


def get_fop_volume(
    fop: Sequence[Sequence[float]], resolution: float, *args: Any, **kwargs: Any
) -> float:
    """Compute field of point's volume.

    Args:
        fop: XYZ coordinates for each atom.
        resolution: resolution in Angstroms.
    """
    return float(len(fop) * (resolution**3))
