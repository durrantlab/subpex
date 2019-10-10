from jdistance import *

if __name__ == "__main__":
    # Get arguments and load json settings file.
    parser = argparse.ArgumentParser(description="Obtain jaccard distance using a reference, a topology file and a MD trajectory file.")
    parser.add_argument("reference", type=str, help="Define the reference PDB file. It is required")
    parser.add_argument("bstate", type=str, help="Define the bstate PDB file. It is required")
    parser.add_argument("settings", type=str, help="Define the json file with the settings. It is required")
    parser.add_argument("--debug", action="store_true", help="To keep all the files to debug any problem")
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
   
    # Define the pocket atoms as all atoms that are at a radius distance of the center point as defined by the user. 
    pocket_reference = reference.select_atoms("point {} {} {} {}".format(str(settings["center"][0]), 
        str(settings["center"][1]), str(settings["center"][2]), str(settings["radius"])))

    # select heavy atoms
    pocket_reference = pocket_reference.select_atoms("not name H")

    # Get the indexes of atoms selected in the reference 
    indexes = pocket_reference.ix

    # Use the indexes to create a selection string to pass selection from reference to ensemble
    selection_pocket = "bynum "+str(indexes[0]) 
    for i in indexes[1:]:
        selection_pocket += " or bynum "+str(i+1)
    #ensemble_pocket = ensemble.select_atoms(selection_pocket)

    # Obtain coordinates for reference atoms and coordinates for alpha carbons.
    reference_coordinates = reference.positions
    reference_alpha = reference.select_atoms("name CA").positions
    reference_fop = get_field_of_points(reference_coordinates, reference_alpha, settings["center"],
                                        settings["resolution"], settings["radius"])


    bb_rmsd = MDAnalysis.analysis.rms.rmsd(reference.select_atoms("backbone").positions, bstate.select_atoms("backbone").positions, 
        center=True, superposition=True)
    pocket_rmsd = MDAnalysis.analysis.rms.rmsd(pocket_reference.positions, bstate.select_atoms(selection_pocket).positions, 
        center=True, superposition=True)
    frame_coordinates = bstate.select_atoms("protein").positions
    frame_calpha = bstate.select_atoms("name CA").positions
    frame_fop = get_field_of_points(frame_coordinates, frame_calpha, settings["center"],
                                      settings["resolution"], settings["radius"])

    pvol = len(frame_fop) * settings['resolution'] ** 3
    jaccard = get_jaccard_distance(reference_fop, frame_fop, settings["resolution"])

    
    print(str(jaccard)+"    "+str(bb_rmsd))

    with open("pvol.txt", "w") as f:
        f.write(str(pvol))


