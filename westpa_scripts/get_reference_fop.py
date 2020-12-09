from jdistance import *

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Obtain FOP for the reference structure")
    parser.add_argument("settings", type=str, help="Define the json file with the settings. It is required")
    args = parser.parse_args()
    try:
        settings = check_input_settings_file(args.settings)
    except IOError:
        print("Could not load the json file with the settings")
        print("make sure the file exists and is correctly formatted")
        raise IOError("Could not load the json file with the settings")

    # Load the reference file defined in the json file
    ref_universe = MDAnalysis.Universe(settings["reference"])
    reference = ref_universe.select_atoms("protein")

    # Define the pocket atoms as all atoms that are at a radius distance of the center point as defined by the user.
    pocket_reference = reference.select_atoms("point {} {} {} {}".format(str(settings["center"][0]),
                                                                         str(settings["center"][1]),
                                                                         str(settings["center"][2]),
                                                                         str(settings["radius"])))

    # obtain all residues in the pocket and create a selection pocket string.
    pocket_residues = []
    for i in pocket_reference:
        residue = str(i.residue)
        residue = residue.split(",")[1][:-1].strip()
        if residue not in pocket_residues:
            pocket_residues.append(residue)
        else:
            pass

    selection_pocket = "resid {} ".format(pocket_residues[0])
    for i in pocket_residues[1:]:
        selection_pocket += " or resid {} ".format(str(i))

    selection_pocket += "and (not name H*)"

    # save the pocket selection string so it can be used in other places
    with open(settings["selection_file"], "w") as f:
        f.write(selection_pocket)

    # use selection pocket string to select the pocket and generate the reference field of points (FOP)
    pocket_reference = reference.select_atoms(selection_pocket)
    reference_coordinates = reference.positions
    reference_alpha = pocket_reference.select_atoms("name CA").positions
    reference_fop = get_field_of_points_dbscan(reference_coordinates, reference_alpha, settings["center"],
                                        settings["resolution"], settings["radius"])

    # Save FOP to a xyz, pdb or pickle file.

    if settings["fop_filetype"] == "xyz":
        points_to_xyz_file(settings["reference_fop"], reference_fop, settings["resolution"], settings["radius"])
    elif settings["fop_filetype"] == "pickle":
        import pickle
        with open(settings["reference_fop"], "wb") as f:
            pickle.dump(reference_fop, f)
    elif settings["fop_filetype"] == "pdb":
        points_to_pdb(settings["reference_fop"], reference_fop)
    else:
        print("Could not save FOP file")
