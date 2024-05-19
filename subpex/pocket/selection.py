import argparse
from collections.abc import Sequence

import MDAnalysis as mda

from ..configs import SubpexConfig
from ..fop.io import write_fop
from .detect import get_fop_pocket


def get_pocket_selection_string(spx_config: SubpexConfig) -> str:
    pocket_selection_mda = spx_config.pocket.selection_mda_str
    pocket_selection_resids = spx_config.pocket.selection_resids
    pocket_selection_dist = spx_config.pocket.selection_dist
    if pocket_selection_mda is not None:
        pocket_selection_str = pocket_selection_mda
    elif pocket_selection_resids is not None:
        pocket_selection_str = "resid " + " or resid ".join(
            [str(i) for i in pocket_selection_resids]
        )
    elif pocket_selection_dist:
        # Define the pocket atoms as all atoms that are at a radius distance of the center point as defined by the user.
        # TODO: Write selection string for this based on distance.
        raise NotImplementedError("Distance pocket selection is not yet implemented")
    else:
        raise RuntimeError(
            "One of MDA, residue IDs, or distance pocket criteria must be specified"
        )
    if spx_config.pocket.selection_append is not None:
        pocket_selection_str += " " + spx_config.pocket.selection_append
    if spx_config.pocket.selection_write is not None:
        with open(spx_config.pocket.selection_write, "w", encoding="utf-8") as f:
            f.write(pocket_selection_str)

    return pocket_selection_str


def get_fop_ref(
    pocket_selection_str: str, spx_config: SubpexConfig
) -> Sequence[Sequence[float]]:
    ref_universe = mda.Universe(spx_config.pocket.ref_pdb_path)
    reference = ref_universe.select_atoms("protein")
    # use selection pocket string to select the pocket and generate the reference field of points (FOP)
    pocket_reference = ref_universe.select_atoms(pocket_selection_str)
    reference_coordinates = reference.positions
    reference_alpha = pocket_reference.select_atoms("name CA").positions
    if spx_config.pocket.center is None:
        raise RuntimeError("The pocket center must be specified.")
    reference_fop = get_fop_pocket(
        protein=reference_coordinates,
        alphas=reference_alpha,
        center=spx_config.pocket.center,
        resolution=spx_config.pocket.resolution,
        radius=spx_config.pocket.radius,
    )
    return reference_fop


def run_pocket_selection(spx_config: SubpexConfig) -> None:
    pocket_selection_str = get_pocket_selection_string(spx_config)
    reference_fop = get_fop_ref(pocket_selection_str, spx_config)
    if spx_config.pocket.ref_fop_write is not None:
        write_fop(reference_fop, spx_config)


def cli_run_detect_pocket_selection():
    parser = argparse.ArgumentParser(
        description="Obtain FOP for the reference structure"
    )
    parser.add_argument(
        "--yaml",
        type=str,
        nargs="+",
        help="Paths to YAML files to use for SubpexContext in decreasing precedence.",
    )
    args = parser.parse_args()

    spx_config = SubpexConfig()
    if args.yaml is not None:
        args.yaml = []
    for yaml_path in reversed(args.yaml):
        spx_config.from_yaml(yaml_path)

    run_pocket_selection(spx_config)
