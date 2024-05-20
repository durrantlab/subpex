from typing import MutableMapping

import MDAnalysis as mda

from subpex.configs import SubpexConfig
from subpex.pocket.detect import get_fop_pocket_convenience, get_pocket_selection


def test_pocket_selection_m7g(
    path_m7g_paths: MutableMapping[str, str], m7g_config: SubpexConfig
) -> None:
    """
    Test the selection of atoms for the binding pocket in a molecular dynamics simulation.

    This function tests the framework's ability to correctly select atoms that define the
    binding pocket based on a specified center, radius, and distance constraint. The selection
    is compared against a known selection string from the configuration.

    Args:
        path_m7g_paths (dict): Dictionary containing paths to the input files required for the test:
            - "ref_pdb": Path to the reference PDB file.
        m7g_config (object): Configuration object for the simulation which includes the parameters
            for defining the pocket (center, radius, and selection distance).

    Raises:
        AssertionError: If the computed pocket selection string does not match the expected value.
    """
    # Load structure
    u = mda.Universe(path_m7g_paths["ref_pdb"])

    pocket_selection = get_pocket_selection(
        u,
        center=m7g_config.pocket.center,  # type: ignore
        radius=m7g_config.pocket.radius,
        distance_constraint=m7g_config.pocket.selection_dist,  # type: ignore
    )
    assert pocket_selection == m7g_config.pocket.selection_str


def test_pocket_fop_gen_m7g(
    path_m7g_paths: MutableMapping[str, str], m7g_config: SubpexConfig
) -> None:
    """
    Test the generation of the Frame of Pocket (FOP) for the binding pocket in a molecular dynamics simulation.

    This function tests the framework's ability to generate a Frame of Pocket (FOP) for the binding
    pocket by detecting the pocket atoms based on the provided configuration. The number of atoms
    and specific coordinates of a known atom are verified against expected values.

    Args:
        path_m7g_paths (dict): Dictionary containing paths to the input files required for the test:
            - "ref_pdb": Path to the reference PDB file.
        m7g_config (object): Configuration object for the simulation which includes settings for
            detecting the pocket.

    Raises:
        AssertionError: If the number of detected pocket atoms or the coordinates of a specific
            atom do not match the expected values.
    """
    u = mda.Universe(path_m7g_paths["ref_pdb"])
    fop_frame = get_fop_pocket_convenience(u.atoms, m7g_config)

    assert len(fop_frame) == 1057
    assert fop_frame[57] == [31.5, 38.5, 31.0]
