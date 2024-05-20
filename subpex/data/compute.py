import argparse
import sys
from collections.abc import MutableMapping, MutableSequence, Sequence

import MDAnalysis as mda
import numpy as np
import pandas as pd
from mda.analysis import align

from ..configs import SubpexConfig
from ..fop.io import read_fop
from .io import initialize_data, write_data


def run_compute_data(
    topo_path: str,
    traj_path: str,
    ref_path: str,
    fop_ref: Sequence[Sequence[float]],
    subpex_config: SubpexConfig,
    selection_align_suffix: str = " and backbone",
    write: bool = True,
    write_dir: str | None = None,
) -> MutableMapping[str, MutableSequence[float | MutableSequence[float]]]:
    """Computes auxillary data during the course of a a simulation.

    Args:
        topo_path: Topology path for the trajectory.
        traj_path: Path to MDAnalysis-supported trajectory to compute progress
            coordinates.
        ref_path: Path to reference structure. This will preferably be a PDB file.
        fop_ref: Reference field of points.
        subpex_config: Context manager for computing progress coordinates.
        selection_align_suffix: Suffix to append to pocket selection for alignment
            purposes.
        write: Write output files.
        write_dir: Directory to write any output files to. Will default to current
            working directory if `None`.
    """
    # Load reference and sample trajectory into universes
    u_ref = mda.Universe(ref_path)
    u_traj = mda.Universe(topo_path, traj_path)

    # Using the selection string to select atoms in the pocket
    if subpex_config.pocket.selection_str is None:
        raise ValueError("pocket_selection_str cannot be None")
    atoms_ref = u_ref.select_atoms("protein")

    # Define the number of times that the progress coordinates and auxiliary
    # info will be calculated.
    num_points = subpex_config.calculated_points
    if num_points == -1:
        num_points = len(u_traj.trajectory)

    # Create a dictionary with all the elements to calculate
    results = initialize_data(subpex_config, data_dir=write_dir)

    # get selection string for alignment
    selection_alignment = subpex_config.pocket.selection_str + selection_align_suffix

    # loop through frames of walker and calculate auxdata
    for frame in np.linspace(0, len(u_traj.trajectory) - 1, num_points, dtype=int):
        u_traj.trajectory[frame]

        # align the frame to the reference
        atoms_frame = u_traj.select_atoms("protein")
        align.alignto(atoms_frame, atoms_ref, select=selection_alignment)

        for aux_data in subpex_config.data.aux:
            if aux_data.active:
                aux_data.compute_frame(
                    atoms_frame=atoms_frame,
                    subpex_config=subpex_config,
                    atoms_ref=atoms_ref,
                )

    write_data(results, subpex_config=subpex_config, data_dir=write_dir)

    return results


def cli_get_data():
    # Get arguments and load json settings file.
    parser = argparse.ArgumentParser(
        description="Compute auxillary data for a trajectory."
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
    args = parser.parse_args()

    subpex_config = SubpexConfig()
    if args.yaml is None:
        args.yaml = []
    for yaml_path in reversed(args.yaml):
        subpex_config.from_yaml(yaml_path)

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

    results = run_compute_data(
        topo_path=args.topo_path,
        traj_path=args.traj_path,
        ref_path=args.ref_path,
        fop_ref=fop_ref,
        subpex_config=subpex_config,
        write_dir=args.write_dir,
    )

    if args.csv is not None:
        df_results = pd.DataFrame(results)
        df_results.to_csv(args.csv)
