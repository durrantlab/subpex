from pcoord import *

def get_pocket_selection(universe, selection_file, fop, distance_constraint):
    with open(settings["selection_file"], "r") as f:
        selection_string = f.readlines()[0]

    pocket_reference = universe.select_atoms(selection_string)

    pocket_list = []
    for atom in pocket_reference:
        if atom.resid in pocket_list:
            continue
        for point in fop:
            distance = calculate_distance_two_points(atom.position, point)
            if distance < distance_constraint:
                if atom.resid not in pocket_list:
                        pocket_list.append(atom.resid)

    selection_pocket = "resid {} ".format(pocket_list[0])
    for i in pocket_list[1:]:
        selection_pocket += " or resid {} ".format(str(i))

    selection_pocket += "and (not name H*)"
    return selection_pocket


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Obtain FOP for the reference structure")
    parser.add_argument("settings", type=str, help="Define the configuration file with the settings. It is required")
    parser.add_argument("dist", type=float, help="Define the distance to the field of points. It is required")
    args = parser.parse_args()
    try:
        settings = check_input_settings_file(args.settings)
    except IOError:
        print("Could not load the json file with the settings")
        print("make sure the file exists and is correctly formatted")
        raise IOError("Could not load the json file with the settings")

    ref_universe = MDAnalysis.Universe(settings["reference"])

    # open file with reference field of points
    if settings["fop_filetype"] == "xyz":
        fop = parse_xyz_fop(settings["reference_fop"])
    elif settings["fop_filetype"] == "pdb":
        reference_fop = parse_pdb_fop(settings["reference_fop"]) #todo make parse_pdb_fop function
    elif settings["fop_filetype"] == "pickle":
        import pickle
        with open(settings["reference_fop"], "rb") as f:
            reference_fop = pickle.load(f)
    else:
        print("could not open reference FOP")
        raise IOError("could not open reference FOP")

    selection_pocket = get_pocket_selection(ref_universe, settings["selection_file"], fop, args.dist)    

    with open(settings["selection_file"], "w") as f:
            f.write(selection_pocket)
