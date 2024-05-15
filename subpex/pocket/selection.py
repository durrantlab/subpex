import argparse
from collections.abc import Sequence

import MDAnalysis as mda

from ..contexts.subpex import SubpexContextManager
from ..fop.io import write_fop
from .detect import get_fop_pocket


def get_pocket_selection_string(subpex_context_manager: SubpexContextManager) -> str:
    pocket_selection_mda = subpex_context_manager.pocket_selection_mda
    pocket_selection_resids = subpex_context_manager.pocket_selection_resids
    pocket_selection_dist = subpex_context_manager.pocket_selection_dist
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
    if subpex_context_manager.pocket_selection_append is not None:
        pocket_selection_str += " " + subpex_context_manager.pocket_selection_append
    if subpex_context_manager.pocket_selection_write is not None:
        with open(
            subpex_context_manager.pocket_selection_write, "w", encoding="utf-8"
        ) as f:
            f.write(pocket_selection_str)

    return pocket_selection_str


def get_fop_ref(
    pocket_selection_str: str, subpex_cm: SubpexContextManager
) -> Sequence[Sequence[float]]:
    ref_universe = mda.Universe(subpex_cm.ref_pdb_path)
    reference = ref_universe.select_atoms("protein")
    # use selection pocket string to select the pocket and generate the reference field of points (FOP)
    pocket_reference = ref_universe.select_atoms(pocket_selection_str)
    reference_coordinates = reference.positions
    reference_alpha = pocket_reference.select_atoms("name CA").positions
    if subpex_cm.pocket_center is None:
        raise RuntimeError("The pocket center must be specified.")
    reference_fop = get_fop_pocket(
        protein=reference_coordinates,
        alphas=reference_alpha,
        center=subpex_cm.pocket_center,
        resolution=subpex_cm.pocket_resolution,
        radius=subpex_cm.pocket_radius,
    )
    return reference_fop


def run_pocket_selection(subpex_cm: SubpexContextManager) -> None:
    pocket_selection_str = get_pocket_selection_string(subpex_cm)
    reference_fop = get_fop_ref(pocket_selection_str, subpex_cm)
    if subpex_cm.ref_fop_write is not None:
        write_fop(reference_fop, subpex_cm)


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

    subpex_cm = SubpexContextManager()
    if args.yaml is None:
        args.yaml = []
    for yaml_path in reversed(args.yaml):
        subpex_cm.from_yaml(yaml_path)

    run_pocket_selection(subpex_cm)
