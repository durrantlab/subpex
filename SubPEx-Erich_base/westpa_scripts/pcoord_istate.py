from jdistance import *

if __name__ == "__main__":
    # Get arguments and load json settings file.
    parser = argparse.ArgumentParser(description="Obtain jaccard distance using a reference, a topology file and a MD trajectory file.")
    parser.add_argument("reference", type=str, help="Define the reference PDB file. It is required")
    parser.add_argument("bstate", type=str, help="Define the bstate PDB file. It is required")
    parser.add_argument("settings", type=str, help="Define the json file with the settings. It is required")
    parser.add_argument("--pvol", action="store_true", help="will save pvol.txt file for auxiliary info in SubPEx run")
    parser.add_argument("--rog", action="store_true", help="will save rog.txt file for auxiliary info in SubPEx run")
    parser.add_argument("--bb_rmsd", action="store_true",
                        help="will save bb_rmsd.txt file for auxiliary info in SubPEx run")
    args = parser.parse_args()

    try:
        with open(args.settings, "r") as f:
            settings = json.load(f)
    except IOError:
        print("Could not load the json file with the settings")
        print("make sure the file exists and is correctly formatted")
        raise IOError("Could not load the json file with the settings")

    # Load reference pdb file and trajectory to then align the trajectory using the reference.
    protein = MDAnalysis.Universe(args.reference)
    reference = protein.select_atoms("protein")
    bstate = MDAnalysis.Universe(args.bstate)

    with open(settings["selection_file"], "r") as f:
        selection_pocket = f.readlines()[0]

    with open(settings["reference_fop"], "rb") as f:
        reference_fop = pickle.load(f)

    # Define the pocket atoms as all atoms that are at a radius distance of the center point as defined by the user.
    pocket_reference = reference.select_atoms(selection_pocket)

    # Obtain coordinates for reference atoms and coordinates for alpha carbons.
    reference_coordinates = reference.positions
    reference_alpha = reference.select_atoms("name CA").positions
    #reference_fop = get_field_of_points(reference_coordinates, reference_alpha, settings["center"],
    #                                    settings["resolution"], settings["radius"])

    pocket_rmsd = MDAnalysis.analysis.rms.rmsd(pocket_reference.positions, bstate.select_atoms(selection_pocket).positions, 
        center=True, superposition=True)
    frame_coordinates = bstate.select_atoms("protein").positions
    frame_calpha = bstate.select_atoms("name CA").positions
    frame_fop = get_field_of_points(frame_coordinates, frame_calpha, settings["center"],
                                      settings["resolution"], settings["radius"])
    jaccard = get_jaccard_distance(reference_fop, frame_fop, settings["resolution"])

    with open("pcoord.txt", "w") as f:
        f.write(str(jaccard)+"    "+str(pocket_rmsd))

    with open("fop.txt", "wb") as f:
        pickle.dump(frame_fop, f)

    if args.bb_rmsd:
        bb_rmsd = MDAnalysis.analysis.rms.rmsd(reference.select_atoms("backbone").positions,
                                               bstate.select_atoms("backbone").positions,
                                               center=True, superposition=True)
        with open("bb_rmsd.txt", "w") as f:
            f.write(str(bb_rmsd))

    if args.pvol:
        pvol = len(frame_fop) * settings['resolution'] ** 3
        with open("pvol.txt", "w") as f:
            f.write(str(pvol))

    if args.rog:
        rog = calculate_pocket_gyration(frame_fop)
        with open("rog.txt", "w") as f:
            f.write(str(rog))
