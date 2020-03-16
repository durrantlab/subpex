from jdistance import *

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Obtain fop for the reference structure")
    parser.add_argument("settings", type=str, help="Define the json file with the settings. It is required")
    args = parser.parse_args()
    try:
        with open(args.settings, "r") as f:
            settings = json.load(f)
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

    with open(settings["selection_file"], "w") as f:
        f.write(selection_pocket)

    reference_fop = get_field_of_points_dbscan(reference_coordinates, reference_alpha, settings["center"],
                                        settings["resolution"], settings["radius"])

    # Save FOP to a xyz file.
    # points_to_xyz_file(settings["reference_fop"], reference_fop, settings["resolution"], settings["radius"])
    # Save fop as a pickle in file given in settings file.
    with open(settings["reference_fop"], "wb") as f:
        pickle.dump(ref_universe, f)
