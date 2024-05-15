from typing import List

import argparse
import sys

import MDAnalysis
import numpy as np
import scipy as sp
from MDAnalysis.analysis import align

from .contexts.pcoord import PCoordContextManager


def get_jaccard_distance(
    reference_fop: List[List[float]], segment_fop: List[List[float]], resolution: float
) -> float:
    """Calculates the Jaccard distance between the points in reference_fop and
    the segment_fop. Uses the distance between points to calculate the intersection.

    Args:
        reference_fop: reference FOP.
        segment_fop: segment FOP.
        resolution: resolution used to create FOP.

    Returns:
        Jaccard distance.
    """

    if len(segment_fop) == 0:
        # sometimes no points are present in the FOP.
        return 1.0

    # Obtaining the trees for both field of points
    reference_tree = sp.spatial.cKDTree(reference_fop)
    segment_tree = sp.spatial.cKDTree(segment_fop)

    # Obtain the points that are at less than resolution/2.5 (aka have the same
    # coordinates)
    clash_indices = reference_tree.query_ball_tree(
        segment_tree, resolution / 2.5, p=2, eps=0
    )

    # Count the points that intersect and convert to float
    intersection = len([x for x in clash_indices if x]) / 1.0

    # Obtain the union of both FOP
    union = float(len(reference_fop) + len(segment_fop) - intersection)

    # Calculate Jaccard distance
    jaccard = 1 - intersection / union

    return jaccard


def run_compute_pcoord(pcoord_context_manager: PCoordContextManager) -> None:
    """Computes progress coordinate for production run.

    Args:
        pcoord_context_manager: Context manager for computing progress coordinates.
    """
    # Load reference pdb file and trajectory to then align the trajectory using
    # the reference.
    ref_universe = MDAnalysis.Universe(settings["reference"])
    reference = ref_universe.select_atoms("protein")
    ensemble = MDAnalysis.Universe(settings["topology"], args.crd_file)

    # open file with reference field of points
    if settings["fop_filetype"] == "xyz":
        reference_fop = parse_xyz_fop(settings["reference_fop"])
    elif settings["fop_filetype"] == "pdb":
        reference_fop = parse_pdb_fop(settings["reference_fop"])
    elif settings["fop_filetype"] == "pickle":
        import pickle

        with open(settings["reference_fop"], "rb") as f:
            reference_fop = pickle.load(f)
    else:
        print("could not open reference FOP")
        raise IOError("could not open reference FOP")

    # Using the selection string to select atoms in the pocket
    pocket_reference = reference.select_atoms(pcoord_context_manager.selection_pocket)

    # Obtain coordinates for reference atoms and coordinates for alpha carbons.
    reference_coordinates = reference.positions

    # Define the number of times that the progress coordinates and auxiliary
    # info will be calculated.
    if settings["calculated_points"] == -1:
        num_points = len(ensemble.trajectory)
    else:
        num_points = settings["calculated_points"]

    # Create a dictionary with all the elements to calculate
    results = {}

    for i in settings["pcoord"]:
        results[i] = []
    for i in settings["auxdata"]:
        results[i] = []
    if "composite" in results.keys():
        if "prmsd" not in results:
            results["prmsd"] = []
        if "bb" not in results:
            results["bb"] = []

    # this section is for the WESTPA analysis tools to work. The first point
    # must be initial point.
    if args.we:
        with open("parent_pcoord.txt", "r") as f:
            initial_pcoords = f.readlines()[-1].split()

        for i, value in enumerate(initial_pcoords):
            results[settings["pcoord"][i]].append(float(value))

        if "jd" in results.keys() and "fops" in settings["auxdata"]:
            with open("parent_fop.txt", "r") as f:
                results["fops"].append(f.readlines())

        if "pvol" in settings["auxdata"]:
            with open("parent_pvol.txt", "r") as f:
                results["pvol"].append(float(f.readlines()[-1]))

        if "rog" in settings["auxdata"]:
            with open("parent_rog.txt", "r") as f:
                results["rog"].append(float(f.readlines()[-1]))

        if "bb" in settings["auxdata"]:
            with open("parent_bb.txt", "r") as f:
                results["bb"].append(float(f.readlines()[-1]))

        if "prmsd" in settings["auxdata"]:
            with open("parent_prmsd.txt", "r") as f:
                results["prmsd"].append(float(f.readlines()[-1]))

        if "jd" in settings["auxdata"]:
            with open("parent_jd.txt", "r") as f:
                results["jd"].append(float(f.readlines()[-1]))

        if "composite" in settings["auxdata"]:
            with open("parent_composite.txt", "r") as f:
                results["composite"].append(float(f.readlines()[-1]))

    # get selction string for alignment
    selection_alignment = selection_pocket + " and backbone"

    # if we are calculating the composite
    if "composite" in results.keys():
        if "sigma" in settings:
            sigma = settings["sigma"]
        else:
            sigma = (
                1
                - len(reference.select_atoms(selection_alignment))
                / len(reference.select_atoms("backbone"))
            ) / 2

    # loop through frames of walker and calculate pcoord and auxdata
    for frame in np.linspace(0, len(ensemble.trajectory) - 1, num_points, dtype=int):
        ensemble.trajectory[frame]
        protein = ensemble.select_atoms("protein")
        # align the frame to the reference
        align.alignto(protein, reference, select=selection_alignment)
        # calculate pocket RMSD if needed
        if "prmsd" in results.keys():
            results["prmsd"].append(
                MDAnalysis.analysis.rms.rmsd(
                    pocket_reference.positions,
                    ensemble.select_atoms(selection_pocket).positions,
                )
            )

        # The next lines calculate the Jaccard distance of the pocket comparing
        # it to the reference
        if (
            "jd" in results.keys()
            or "fops" in results.keys()
            or "pvol" in results.keys()
            or "rog" in results.keys()
        ):
            frame_coordinates = ensemble.select_atoms("protein").positions
            pocket_calpha = ensemble.select_atoms(
                selection_pocket + " and name CA*"
            ).positions
            frame_fop = get_fop_pocket(
                frame_coordinates,
                pocket_calpha,
                settings["center"],
                settings["resolution"],
                settings["radius"],
            )
            results["jd"].append(
                get_jaccard_distance(reference_fop, frame_fop, settings["resolution"])
            )

        if "fops" in results.keys():
            results["fops"].append(frame_fop)

        if "pvol" in results.keys():
            results["pvol"].append(len(frame_fop) * (settings["resolution"] ** 3))

        if "rog" in results.keys():
            results["rog"].append(calculate_pocket_gyration(frame_fop))

        if "bb" in results.keys():
            align.alignto(protein, reference, select="backbone")
            results["bb"].append(
                MDAnalysis.analysis.rms.rmsd(
                    reference.select_atoms("backbone").positions,
                    protein.select_atoms("backbone").positions,
                )
            )

        if "composite" in results.keys():
            results["composite"].append(
                results["prmsd"][-1] + (sigma * results["bb"][-1])
            )

    # writing in text files the progress coordinates and the required auxiliary
    # data if needed.
    with open("pcoord.txt", "w") as f:
        for i in range(len(results[settings["pcoord"][0]])):
            line = ""
            for pcoord in settings["pcoord"]:
                line += "{:.4f}    ".format(results[pcoord][i])
            f.write(line + "\n")

    # save fop in file so it can be piped to h5 file
    if "fops" in results.keys():
        if settings["fop_filetype"] == "xyz":
            points_to_xyz(
                "fop.txt", results["fops"], settings["resolution"], settings["radius"]
            )
        elif settings["fop_filetype"] == "pdb":
            points_to_pdb("fop", results["fops"])
        # not sure pickles work
        elif settings["fop_filetype"] == "pickle":
            with open("fop.txt", "wb") as f:
                pickle.dump(frame_fop, f)

    if "pvol" in settings["auxdata"]:
        with open("pvol.txt", "w") as f:
            for i in results["pvol"]:
                f.write(str(i) + "\n")

    if "rog" in settings["auxdata"]:
        with open("rog.txt", "w") as f:
            for i in results["rog"]:
                f.write(str(i) + "\n")

    if "bb" in settings["auxdata"]:
        with open("bb.txt", "w") as f:
            for i in results["bb"]:
                f.write(str(i) + "\n")

    if "prmsd" in settings["auxdata"]:
        with open("prmsd.txt", "w") as f:
            for i in results["prmsd"]:
                f.write(str(i) + "\n")

    if "jd" in settings["auxdata"]:
        with open("jd.txt", "w") as f:
            for i in results["jd"]:
                f.write(str(i) + "\n")

    if "composite" in settings["auxdata"]:
        with open("composite.txt", "w") as f:
            for i in results["composite"]:
                f.write(str(i) + "\n")

    if args.csv is not None:
        import pandas as pd

        results_pd = pd.DataFrame(results)
        results_pd.to_csv(args.csv)


def cli_get_pcoord():
    # Get arguments and load json settings file.
    parser = argparse.ArgumentParser(
        description="Obtain jaccard distance using a reference, a topology file and a MD trajectory file."
    )
    parser.add_argument(
        "crd_file", type=str, help="Define the coordinate file. It is required"
    )
    parser.add_argument(
        "--yaml",
        type=str,
        nargs="+",
        help="Paths to YAML files to use for PCoordContext in decreasing precedence.",
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
        "calculation in a weighted ensemble simulation",
    )
    args = parser.parse_args()

    pcoord_context_manager = PCoordContextManager()
    if args.yaml is None:
        args.yaml = []
    for yaml_path in reversed(args.yaml):
        pcoord_context_manager.from_yaml(yaml_path)

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

    run_compute_pcoord(pcoord_context_manager)
