import argparse
import sys
from collections.abc import MutableMapping, Sequence

import MDAnalysis as mda
import numpy as np
import pandas as pd
from mda.analysis import align

from ..contexts.subpex import SubpexContextManager
from ..fop.io import read_fop
from ..fop.prop import get_fop_volume
from ..pocket.detect import get_fop_pocket
from ..pocket.props import get_pocket_rog
from .desc import get_rmsd
from .io import initialize_data, write_data
from .jaccard import get_jaccard_distance


def compute_frame_data(
    fop_ref: Sequence[Sequence[float]],
    atoms_frame: mda.AtomGroup,
    atoms_ref: mda.AtomGroup,
    atoms_ref_pocket: mda.AtomGroup,
    subpex_cm: SubpexContextManager,
    composite_sigma: float | None = None,
) -> MutableMapping[str, float]:
    results_frame: MutableMapping[str, float] = {}
    data_keys = subpex_cm.data_keys

    if subpex_cm.key_pocket_rmsd in data_keys:
        results_frame[subpex_cm.key_pocket_rmsd] = get_rmsd(
            atoms_ref_pocket.positions,
            atoms_frame.select_atoms(subpex_cm.pocket_selection_str).positions,
        )

    # The next lines calculate the Jaccard distance of the pocket comparing
    # it to the reference
    if (
        "jd" in data_keys
        or "fops" in data_keys
        or "p_vol" in data_keys
        or "p_rog" in data_keys
    ):
        frame_coordinates = atoms_frame.select_atoms("protein").positions
        if subpex_cm.pocket_selection_str is None:
            raise ValueError("Pocket selection string cannot be `None`")
        pocket_calpha = atoms_frame.select_atoms(
            subpex_cm.pocket_selection_str + " and name CA*"
        ).positions
        if subpex_cm.pocket_center is None:
            raise ValueError("`pocket_center` cannot be None")
        fop_frame = get_fop_pocket(
            frame_coordinates,
            pocket_calpha,
            subpex_cm.pocket_center,
            subpex_cm.pocket_resolution,
            subpex_cm.pocket_radius,
        )

    if subpex_cm.key_pocket_jd in data_keys:
        results_frame[subpex_cm.key_pocket_jd] = get_jaccard_distance(
            fop_ref, fop_frame, subpex_cm
        )

    if subpex_cm.key_pocket_fop in data_keys:
        results_frame[subpex_cm.key_pocket_fop] = fop_frame

    if subpex_cm.key_pocket_vol in data_keys:
        results_frame[subpex_cm.key_pocket_vol] = get_fop_volume(
            fop_frame, subpex_cm.pocket_resolution
        )

    if subpex_cm.key_pocket_rog in data_keys:
        results_frame[subpex_cm.key_pocket_rog] = get_pocket_rog(fop_frame)

    if subpex_cm.key_bb_rmsd in data_keys:
        align.alignto(atoms_frame, atoms_ref, select="backbone")
        results_frame[subpex_cm.key_bb_rmsd] = get_rmsd(
            atoms_ref.select_atoms("backbone").positions,
            atoms_frame.select_atoms("backbone").positions,
        )

    if subpex_cm.key_composite in data_keys:
        if subpex_cm.key_pocket_rmsd not in data_keys:
            raise RuntimeError("pocket_rmsd must be active for composite")
        if subpex_cm.key_pocket_rmsd not in data_keys:
            raise RuntimeError("backbone_rmsd must be active for composite")
        if composite_sigma is None:
            raise RuntimeError("composite_sigma cannot be None for composite")
        results_frame[subpex_cm.key_composite] = results_frame[
            subpex_cm.key_pocket_rmsd
        ] + (composite_sigma * results_frame[subpex_cm.key_bb_rmsd])

    return results_frame


def run_compute_pcoord(
    topo_path: str,
    traj_path: str,
    ref_path: str,
    fop_ref: Sequence[Sequence[float]],
    subpex_cm: SubpexContextManager,
    selection_align_suffix: str = " and backbone",
    write: bool = True,
    write_dir: str | None = None,
) -> MutableMapping[str, MutableSequence[float | MutableSequence[float]]]:
    """Computes progress coordinate for production run.

    Args:
        topo_path: Topology path for the trajectory.
        traj_path: Path to MDAnalysis-supported trajectory to compute progress
            coordinates.
        ref_path: Path to reference structure. This will preferably be a PDB file.
        fop_ref: Reference field of points.
        subpex_cm: Context manager for computing progress coordinates.
        selection_align_suffix: Suffix to append to pocket selection for alignment
            purposes.
        write: Write output files.
        write_dir: Directory to write any output files to. Will default to current
            working directory if `None`.
    """
    # Load reference and sample trajectory into universes
    u_ref = mda.Universe(ref_path)
    u_ensemble = mda.Universe(topo_path, traj_path)

    # Using the selection string to select atoms in the pocket
    reference = u_ref.select_atoms("protein")
    pocket_reference = reference.select_atoms(subpex_cm.pocket_selection_str)

    # Define the number of times that the progress coordinates and auxiliary
    # info will be calculated.
    num_points = subpex_cm.calculated_points
    if num_points == -1:
        num_points = len(u_ensemble.trajectory)

    # Create a dictionary with all the elements to calculate
    results = initialize_data(subpex_cm, data_dir=write_dir)

    # get selection string for alignment
    if subpex_cm.pocket_selection_str is None:
        raise ValueError("pocket_selection_str cannot be None")
    selection_alignment = subpex_cm.pocket_selection_str + selection_align_suffix

    # if we are calculating the composite
    if "composite" in results.keys():
        if subpex_cm.composite_sigma is None:
            composite_sigma = (
                1
                - len(reference.select_atoms(selection_alignment))
                / len(reference.select_atoms("backbone"))
            ) / 2
        else:
            composite_sigma = subpex_cm.composite_sigma

    # loop through frames of walker and calculate pcoord and auxdata
    for frame in np.linspace(0, len(u_ensemble.trajectory) - 1, num_points, dtype=int):
        u_ensemble.trajectory[frame]

        # align the frame to the reference
        protein = u_ensemble.select_atoms("protein")
        align.alignto(protein, reference, select=selection_alignment)

        results_frame = compute_frame_data(
            fop_ref=fop_ref,
            atoms_frame=protein,
            atoms_ref=reference,
            atoms_ref_pocket=pocket_reference,
            subpex_cm=subpex_cm,
            composite_sigma=composite_sigma,
        )
        for key, value in results_frame.items():
            results[key].append(value)

    write_data(results, subpex_cm=subpex_cm, data_dir=write_dir)

    return results


def cli_get_pcoord():
    # Get arguments and load json settings file.
    parser = argparse.ArgumentParser(
        description="Obtain jaccard distance using a reference, a topology file and a MD trajectory file."
    )
    parser.add_argument("topo_path", type=str, help="Path to topology file.")
    parser.add_argument(
        "ref_path", type=str, help="Path to reference structure (preferably PDB)."
    )
    parser.add_argument(
        "traj_path", type=str, help="Path to trajectory file supported by mda."
    )
    parser.add_argument("fop_ref_path", type=str, help="Path to reference FOP file.")
    parser.add_argument(
        "--write_dir",
        type=str,
        help="Directory to write, and search for, any output files.",
    )
    parser.add_argument(
        "--yaml",
        type=str,
        nargs="+",
        help="Paths to YAML files to use for subpex context in decreasing precedence.",
    )
    parser.add_argument(
        "--csv",
        type=str,
        help="will save results in an csv file with the provided filename",
    )
    parser.add_argument(
        "--we",
        action="store_true",
        help="Use this when you are using this for a progress coordinate "
        "calculation in a weighted u_ensemble simulation",
    )
    args = parser.parse_args()

    subpex_cm = SubpexContextManager()
    if args.yaml is None:
        args.yaml = []
    for yaml_path in reversed(args.yaml):
        subpex_cm.from_yaml(yaml_path)

    # Checking the the trajectory doesn't have a rattle error.
    try:
        with open("seg.log", "r", encoding="utf-8") as f:
            line = f.readline()
            while line:
                if "RATTLE" in line:
                    print("There was an error with the simulation")
                    sys.exit("There was an error with the simulation")
                else:
                    line = f.readline()
    except Exception:
        print("Could not check log file. There could be a problem with the simulation")

    fop_ref = read_fop(args.fop_ref_path)

    results = run_compute_pcoord(
        topo_path=args.topo_path,
        traj_path=args.traj_path,
        ref_path=args.ref_path,
        fop_ref=fop_ref,
        subpex_cm=subpex_cm,
        write_dir=args.write_dir,
    )

    if args.csv is not None:
        df_results = pd.DataFrame(results)
        df_results.to_csv(args.csv)
