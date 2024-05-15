from pcoord import *

# todo refactor this script
if __name__ == "__main__":
    # Get arguments and load json settings file.
    parser = argparse.ArgumentParser(
        description="Obtain jaccard distance using a reference, a topology file and a MD trajectory file."
    )
    parser.add_argument(
        "istate", type=str, help="Define the istate or istate PDB file. It is required"
    )
    parser.add_argument(
        "settings", type=str, help="Define the settings file. It is required"
    )
    args = parser.parse_args()

    # Function that obtains and checks for settings in the settings json file.
    print(args.settings)
    settings = check_input_settings_file(args.settings)

    # Load reference pdb file and trajectory to then align the trajectory using the reference.
    protein = MDAnalysis.Universe(settings["reference"])
    reference = protein.select_atoms("protein")

    istate = MDAnalysis.Universe(args.istate)

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

    # open file with selection string
    with open(settings["selection_file"], "r") as f:
        selection_pocket = f.readlines()[0]

    selection_alignment = selection_pocket + " and backbone"

    # if we are calculating composite progress coordiante (or auxiliary data) check that we have sigma
    if "composite" in results.keys():
        if "sigma" in settings:
            sigma = settings["sigma"]
        else:
            sigma = (
                1
                - len(reference.select_atoms(selection_alignment))
                / len(reference.select_atoms("backbone"))
            ) / 2

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

    align.alignto(istate, reference, select=selection_alignment)

    # Using the selection string to select atoms in the pocket
    pocket_reference = reference.select_atoms(selection_pocket)

    # calculate pocket rmsd

    if "prmsd" in results.keys():
        results["prmsd"].append(
            MDAnalysis.analysis.rms.rmsd(
                pocket_reference.positions,
                istate.select_atoms(selection_pocket).positions,
            )
        )

    if (
        "jd" in results.keys()
        or "fops" in results.keys()
        or "pvol" in results.keys()
        or "rog" in results.keys()
    ):
        frame_coordinates = istate.select_atoms("protein").positions
        pocket_calpha = istate.select_atoms(
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
        align.alignto(istate, reference, select="backbone")
        results["bb"].append(
            MDAnalysis.analysis.rms.rmsd(
                reference.select_atoms("backbone").positions,
                istate.select_atoms("backbone").positions,
            )
        )

    if "composite" in results.keys():
        results["composite"].append(results["prmsd"][-1] + (sigma * results["bb"][-1]))

    # writing in text files the progress coordinates and the required auxiliary data if needed.
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
