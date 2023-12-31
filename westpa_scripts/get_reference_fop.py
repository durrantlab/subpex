from pcoord import *
from typing import List


def get_pocket_selection(
    universe,
    center: List[float],
    radius: float,
    selection_file: str,
    distance_constraint: float = 6.7,
) -> str:
    """Takes the protein and center of the pocket and gets the initial
    selection of the pockets based on the user-specified radius. If water
    molecules are present, it uses them as a first approximation of the surface
    residues.

    Args:
        universe (MDAnalysis Universe): MDA universe with protein 
        center (list): x, y, z coordiantes of the pocket center
        radius (float): radius of the pocket to be considered
        selection_file (str): filename for the selection pocket
        distance_constraint (float, optional): Constraint to use for proximity
            to calculate surface residues. Defaults to 6.7.

    Returns:
        selection_pocket (str): selection string of the pocket as used in
            MDAnalysis
    """

    water_pocket = universe.select_atoms(
        f"point {center[0]} {center[1]} {center[2]} {radius} and resname WAT"
    )
    residues_pocket = universe.select_atoms(
        f"protein and point {center[0]} {center[1]} {center[2]} {radius}"
    )

    if len(water_pocket.residues) == 0:
        print(
            "There are no water molecules to obtain distances, some burried residues may end up in the pocket selection"
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
        # we want to make sure that the residues we are selecting are close to the surface, we will approximate the surface accesible residues
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

        # creating selection pocket string comaptible with
        selection_pocket = f"resid {list_close_water[0]} "
        for i in list_close_water[1:]:
            selection_pocket += f" or resid {str(i)} "

            # save the pocket selection string so it can be used in other places
    selection_pocket += "and (not name H*)"

    with open(selection_file, "w") as f:
        f.write(selection_pocket)

    return selection_pocket


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Obtain FOP for the reference structure"
    )
    parser.add_argument(
        "settings",
        type=str,
        help="Define the yaml file with the settings. Most likely should be west.cfg. It is required",
    )
    parser.add_argument(
        "-mda",
        "--selection_mda",
        type=str,
        help="Define selection for the pocket atoms. It has to a selection string compatible with MDAnalysis",
    )
    parser.add_argument(
        "-resid",
        "--selection_resid",
        type=str,
        help="Define selection for the pocket atoms. It has to a selection list of residue numebers (e.g. 1 2 3 4 5 6)",
    )
    parser.add_argument(
        "-dist",
        "--distance_constraint",
        type=float,
        help="Define the distance cosntraint between water and residues center of geometry for pocket selection. Default 6.7",
    )
    args = parser.parse_args()

    try:
        settings = check_input_settings_file(args.settings)
    except IOError as e:
        print("Could not load the json file with the settings")
        print("make sure the file exists and is correctly formatted")
        raise IOError("Could not load the json file with the settings") from e

    # Load the reference file defined in the json file
    ref_universe = MDAnalysis.Universe(settings["reference"])
    reference = ref_universe.select_atoms("protein")

    # Define the pocket by the users selection
    if args.selection_mda:
        try:
            pocket_reference = reference.select_atoms(args.selection_mda)
            selection_pocket = args.selection_mda
        except ValueError:
            sys.exit("There is an error on your selection string")
        with open(settings["selection_file"], "w") as f:
            f.write(args.selection_mda)

    elif args.selection_resid:
        # Define the pocket by the user list of residue numbers
        try:
            [int(x) for x in args.selection_resid.split()]
        except ValueError:
            sys.exit("Make sure the list only contains residue numbers")

        selection_pocket = f"resid {args.selection_resid.split()[0]} "
        for i in args.selection_resid.split()[1:]:
            selection_pocket += f" or resid {str(i)} "
        selection_pocket += "and (not name H*)"
        with open(settings["selection_file"], "w") as f:
            f.write(selection_pocket)

    elif args.distance_constraint:
        # Define the pocket atoms as all atoms that are at a radius distance of the center point as defined by the user.
        selection_pocket = get_pocket_selection(
            ref_universe,
            settings["center"],
            settings["radius"],
            settings["selection_file"],
            args.distance_constraint,
        )
    else:
        selection_pocket = get_pocket_selection(
            ref_universe,
            settings["center"],
            settings["radius"],
            settings["selection_file"],
        )

    # use selection pocket string to select the pocket and generate the reference field of points (FOP)
    pocket_reference = ref_universe.select_atoms(selection_pocket)
    reference_coordinates = reference.positions
    reference_alpha = pocket_reference.select_atoms("name CA").positions
    reference_fop = get_field_of_points_dbscan(
        reference_coordinates,
        reference_alpha,
        settings["center"],
        settings["resolution"],
        settings["radius"],
    )

    # Save FOP to a xyz, pdb or pickle file.

    if settings["fop_filetype"] == "xyz":
        points_to_xyz_file(
            settings["reference_fop"],
            reference_fop,
            settings["resolution"],
            settings["radius"],
        )
    elif settings["fop_filetype"] == "pickle":
        import pickle

        with open(settings["reference_fop"], "wb") as f:
            pickle.dump(reference_fop, f)
    elif settings["fop_filetype"] == "pdb":
        points_to_pdb(settings["reference_fop"], reference_fop)
    else:
        print("Could not save FOP file")

    print("")
    print("Files create successfully:")
    print("   " + settings["selection_file"])
    print("   " + settings["reference_fop"])

