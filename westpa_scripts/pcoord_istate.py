from jdistance import *

#todo refactor this script
if __name__ == "__main__":
    # Get arguments and load json settings file.
    parser = argparse.ArgumentParser(description="Obtain jaccard distance using a reference, a topology file and a MD trajectory file.")
    parser.add_argument("istate", type=str, help="Define the istate or istate PDB file. It is required")
    parser.add_argument("settings", type=str, help="Define the settings file. It is required")
    args = parser.parse_args()

    # Function that obtains and checks for settings in the settings json file.
    print(args.settings)
    settings = check_input_settings_file(args.settings)

    # Load reference pdb file and trajectory to then align the trajectory using the reference.
    protein = MDAnalysis.Universe(settings["reference"])
    reference = protein.select_atoms("protein")

    istate = MDAnalysis.Universe(args.istate)

    # open file with selection string
    with open(settings["selection_file"], "r") as f:
        selection_pocket = f.readlines()[0]

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
    pocket_reference = reference.select_atoms(selection_pocket)

    # Obtain coordinates for reference atoms and coordinates for alpha carbons.
    reference_coordinates = reference.positions

    # obtain pocket for istate or bstate and the CA
    istate_pocket = istate.select_atoms(selection_pocket)
    frame_coordinates = istate.select_atoms("protein").positions
    frame_calpha = istate_pocket.select_atoms("name CA").positions

    # calculate pocket rmsd
    pocket_rmsd = MDAnalysis.analysis.rms.rmsd(pocket_reference.positions, istate_pocket.positions,
                                               center=True, superposition=True)

    frame_fop = get_field_of_points_dbscan(frame_coordinates, frame_calpha, settings["center"],
                                      settings["resolution"], settings["radius"])
    jaccard = get_jaccard_distance(reference_fop, frame_fop, settings["resolution"])

    with open("pcoord.txt", "w") as f:
        f.write(str(jaccard)+"    "+str(pocket_rmsd))

    # save fop in file so it can be piped to h5 file
    if settings["fop_filetype"] == "xyz":
        points_to_xyz_file("fop.txt", frame_fop, settings["resolution"], settings["radius"])
    elif settings["fop_filetype"] == "pdb":
        points_to_pdb("fop", results["fops"])
    elif settings["fop_filetype"] == "pickle":
        with open("fop.txt", "wb") as f:
            pickle.dump(frame_fop, f)


    if "bb_rmsd" in settings["auxdata"]:
        bb_rmsd = MDAnalysis.analysis.rms.rmsd(reference.select_atoms("backbone").positions,
                                               istate.select_atoms("backbone").positions,
                                               center=True, superposition=True)
        with open("bb_rmsd.txt", "w") as f:
            f.write(str(bb_rmsd))

    if "pvol" in settings["auxdata"]:
        pvol = len(frame_fop) * settings['resolution'] ** 3
        with open("pvol.txt", "w") as f:
            f.write(str(pvol))

    if "rog" in settings["auxdata"]:
        rog = calculate_pocket_gyration(frame_fop)
        with open("rog.txt", "w") as f:
            f.write(str(rog))
