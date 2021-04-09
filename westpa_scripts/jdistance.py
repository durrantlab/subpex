import argparse
import MDAnalysis
import numpy as np
import scipy as sp
from MDAnalysis.analysis import rms, align
from sklearn.cluster import DBSCAN
import sys
import logging


def check_input_settings_file(filename):
    """
    check_input_settings_file is a function that loads the json file containing the settings for SubPEx and checks that 
    all the parameters are available or sets defaults if it can be done.
    :param filename: (str) filename of the settings json file. 
    :return: (dict) settings for running SubPEx and some of the analysis tools  
    """

    import sys
    import glob

    import yaml
    try:
        with open(filename, "r") as f:
            settings = yaml.safe_load(f.read())
        settings = settings["subpex"]
    except IOError:
        print("Could not load the configuration file with the settings")
        print("make sure the file exists and is correctly formatted")
        raise IOError("Could not load the json file with the settings")

    # Checking that the coordinates of the center of the pocket are in settings.
    if "center" in settings:
        if type(settings["center"]) == list and len(settings["center"]) == 3:
            if type(settings["center"][0]) == float and type(settings["center"][0]) == float and type(
                    settings["center"][0]) == float:
                pass
            else:
                logging.critical(
                    "There is an error with the center parameter in the settings file. Check if each element is a float")
                sys.exit("Error with setting center")
        else:
            logging.critical(
                "There is an error with the center parameter in the settings file. Check if it is a list with the XYZ coordinates")
            sys.exit("Error with settings center")

    # Checking that we specify which progress coordinates we will print in the pcoord file
    if "pcoord" in settings:
        available_pcoords = ["jd", "prmsd", "bb_rmsd", "composite", "pvol", "rog_pocket"]
        if len(settings["pcoord"]) < 1:
            logging.critical("There is an error with the progress coordinate (pcoord) setting. Need to have jd, bb_rmsd and/or prsmd")
            sys.exit("Error with setting pcoord")
        else:
            for i in settings["pcoord"]:
                if i not in available_pcoords:
                    logging.critical("There is an error with the progress coordinate (pcoord) setting. Need to have jd, bb_rmsd and/or prsmd")
                    sys.exit("Error with setting pcoord")
    else:
        logging.critical("There is an error with the progress coordinate (pcoord) setting. Need to have jd, bb_rmsd and/or prsmd")
        sys.exit("Error with setting pcoord")

    # Checking radius is in settings
    if "radius" in settings and type(settings["radius"]) == float:
        pass
    else:
        # The default value of 10.9 Angstrom comes from the most common volume for a pocket accordintg to https://doi.org/10.1002/pro.5560070905 and approximating the pocket as a spherical pocket.
        settings["radius"] = 10.90
        logging.warning(
            "Using a default radius value of 10.90 Angstrom. Please check it will encompases all of your pocket.")

    # Cheking resolution setting
    if "resolution" in settings and type(settings["resolution"]) == float:
        pass
    else:
        settings["resolution"] = 0.5
        logging.warning("Using a default resolution value of 0.5 Angstrom.")

    # Cheking calculated_points setting
    if "calculated_points" in settings and type(settings["calculated_points"]) == int:
        pass
    else:
        settings["calculated_points"] = -1
        logging.warning(
            "Using a default calculated_points value of -1. This will calculate progress coordinates for all the frames in the trajectory file")

    # Checking west_home setting
    if "west_home" not in settings or len(glob.glob(settings["west_home"])) != 1:
        logging.critical("There is an error with the west home directory")
        sys.exit("There is an error with the west home directory")
    elif len(glob.glob(settings["west_home"] + "/traj_segs")) != 1:
        logging.warning(
            "Could not locate the traj_segs directory. Please be sure it exists and the trajectories are there before running the clustering script. Note: This should not exist if you have not run SubPEx")
    else:
        pass

    # checking reference and topology files
    check_file_exists(settings, "reference")
    check_file_exists(settings, "topology")

    # checking that the selection_file exists, it is different for get_reference_fop.py script
    if __file__ == "jdistance.py" or __file__ == "clustering.py":
        check_file_exists(settings, "selection_file")
    else:
        if "selection_file" not in settings:
            settings["selection_file"] = settings["west_home"] + "/selection_string.txt"
            logging.warning("Selection file was not determined used the default of {}/selection_string.txt".format(
                settings["west_home"]))

    # making sure the reference_fop exists and the format is in the settings file
    accepted_fop_files = ["xyz", "pdb"]
    if "fop_filetype" not in settings and settings["fop_filetype"] not in accepted_fop_files:
        settings["fop_filetype"] = "xyz"

    if __file__ == "jdistance.py":
        check_file_exists(settings, "reference_fop")
    else:
        if "reference_fop" not in settings:
            settings["reference_fop"] = settings["west_home"] + "/reference_fop.{}".format(settings["fop_filetype"])
            logging.warning("Selection file was not determined used the default of {}/reference_fop.{}".format(
                settings["west_home"]), settings["fop_filetype"])

    logging.info("The parameters used for {} are {}".format(__file__, settings))
    return settings


def check_file_exists(settings, keyword):
    """
    A simple function that checks that the keyword is in the dictionary settings and that the file associated with the
    keyword exists.
    :param settings: dic dictionary with settings for subpex
    :param keyword: str keyword of the setting to check the file exists.
    :return:
    """
    import sys
    import glob

    if keyword not in settings or len(glob.glob(settings[keyword])) != 1:
        logging.critical("There is an error with the {} parameter, could not locate the file".format(keyword))
        sys.exit("There is an error with the {} file".format(keyword))
    else:
        pass


def points_to_pdb(filename, coordinates):
    """
    points_to_pdb will write a pdb file full of C-alphas to be able to load the file in a visualisation software and
    represent the field of points.

    :param filename: (str) name of the pdb file to create.
    :param coordinates: (list of lists) XYZ coordinates of the field of points to write into pdb file.
    """
    with open(filename, "w") as f: #todo modify to be able to do multiframe pdb
        # f.write(header)
        atom_number = 1
        for i in coordinates:
            text = "ATOM" + "{:>7}".format(atom_number)
            text += "{:^6}".format("CA")
            text += "{:>3}".format("ALA")
            text += "{:>6}".format("A")
            text += "{:>12}".format(i[0])
            text += "{:>8}".format(i[1])
            text += "{:>8}".format(i[2])
            text += "{:>6}".format("1.0")
            text += "{:>6}".format("1.0")
            text += "\n"
            f.write(text)
            atom_number += 1
        f.write("\n")


def parse_xyz_fop(filename):
    """

    :param filename:
    :return:
    """
    # open reference fop xyz file
    with open(filename, "r") as f:
        text = f.readlines()

    # parse xyz file
    fop = []
    for i in text[2:]:
        line = i.split()
        if len(line) == 4:
            point = [float(line[1]), float(line[2]), float(line[3])]
            fop.append(point)
        else:
            pass

    return fop


def points_to_xyz_file(filename, coordinates, resolution, radius):
    """
    points_to_xyz_file takes the coordinates and the resolution and coordinates to write an xyz file that can be read by
    a visualization software.

    :param filename: (str) name for the xyz file to be created.
    :param coordinates: (list of lists) contains XYZ coordinates for each atom.
    :param float resolution: Resolution in Angstroms.
    :param float radius: radius in Angstrom for FOP.
    """
    # calculate the volume using the resolution
    volume = len(coordinates) * (resolution ** 3) #todo modify to be able to do multiframe xyz

    # Writing the xyz file.
    with open(filename, "w") as f:
        f.write(str(len(coordinates)) + "\n")
        f.write("Generated in SubPEx. Res: {res}    Rad: {dim} Ang    Volume: {vol} Ang^3 \n".format(res=resolution,
                                                                                                     dim=radius,
                                                                                                     vol=volume))
        # Adding each coordinate as a C-alpha for visualization purposes
        for i in coordinates:
            f.write("CA    {}    {}    {} \n".format(i[0], i[1], i[2]))


def calculate_distance_two_points(point1, point2):
    """
    This function calculates the distance between two points in a 3D coordinate system.

    :param point1: (list) coordinates of first point. [x, y, z]
    :param point2: (list) coordinates of second point. [x, y, z]
    :return distance: (float) distance between point1 and point2.
    """
    distance = 0
    for i in range(len(point1)):
        distance += (point1[i] - point2[i]) ** 2

    return np.sqrt(distance)


def field_of_points(center, resolution, radius):
    """
    Function that snaps the provided center to the a center congruent to the resolution and radius given by the
    user. Then it generates a spherical field of points (FOP) with radius = radius/2 with points at with a distance
    between them equal to the resolution.

    :param list center: provided center coordinates. [X, Y, Z]
    :param float resolution: Distance between points which give resolution of the pocket volume calculator method.
    :param float radius: Defines the radius of the sphere of the FOP.
    :return list of lists : a spherical FOP
    """
    # This part of the code will get a new center, a snapped center.
    snapped_x = np.around(center[0] / resolution) * resolution
    snapped_y = np.around(center[1] / resolution) * resolution
    snapped_z = np.around(center[2] / resolution) * resolution
    center = np.around(np.array([snapped_x, snapped_y, snapped_z], float), 3)

    # Getting the field of points around the center.
    xs = np.round(np.arange(center[0] - radius, center[0] + radius + resolution, resolution), 3)
    ys = np.round(np.arange(center[1] - radius, center[1] + radius + resolution, resolution), 3)
    zs = np.round(np.arange(center[2] - radius, center[2] + radius + resolution, resolution), 3)
    field_of_points = np.vstack(np.meshgrid(xs, ys, zs)).reshape(3, -1).T

    # making the field of points a sphere of radius = radius/2.
    field_of_points = [x for x in field_of_points if calculate_distance_two_points(x, center) < radius]

    return field_of_points, center


def get_trimmed_fop(field_of_points, atoms_points):
    """
    Removes points that clash with atoms within the FOP.

    :param field_of_points: List of list containing all the points for hte FOP.
    :param atoms_points: List of list containing the coordinates of atoms within the FOP.
    :return: Returns a list of list that has the trimmed FOP.
    """
    # Generating cKDTrees to obtain distances
    atoms_tree = sp.spatial.cKDTree(atoms_points)
    fop_tree = sp.spatial.cKDTree(field_of_points)

    # get indices of fop_tree that are close to atoms_tree
    # We use the Van der Waals radius of water and that of hydrogen atom
    vdw_radius = 2.6
    clash_indices = atoms_tree.query_ball_tree(fop_tree, vdw_radius, p=2, eps=0)
    for index in clash_indices:
        for item in index:
            field_of_points[item] = None

    # Creating the FOP that does not contain points that clash with the atoms provided to the function.
    trimmed_fop = []
    for i in field_of_points:
        if type(i) == np.ndarray:
            trimmed_fop.append(i)
        else:
            pass

    return trimmed_fop


def get_field_of_points_dbscan(protein, alphas, center, resolution, radius):
    """
    This function will take the coordinates of the protein and the coordinates of alpha carbons calculate the field of
    points (FOP).

    :param protein: (list of list) coordinates of protein atoms.
    :param alphas: (list of list) coordinates of all alpha carbons
    :param center: (list) center of the pocket (user defined).
    :param resolution: (float) resolution for the FOP (user defined).
    :param radius: (float) radius of sphere for FOP (user defined).
    :return: list of list with the pocket shape as a FOP.
    """
    # get starting FOP and snapped center
    fop, center = field_of_points(center, resolution, radius)
    # trim protein atoms to those within the FOP sphere.
    trimmed_coords_protein = [x for x in protein if calculate_distance_two_points(x, center) < radius]
    # remove points in FOP that have steric clashes with protein.
    trimmed_fop = get_trimmed_fop(fop, trimmed_coords_protein)
    # remove points outside the convex hull created by the alpha carbons in the pocket.
    pocket = remove_convex_fop(trimmed_fop, alphas)
    pocket = cluster_dbscan(pocket)

    return pocket


def cluster_dbscan(fop):
    """
    This function will take the coordinates of the field of points and perform a clustering.py algorithm to trim the fop
    to the points in the same cluster that the center of the pocket.

    :param protein: (list of list) coordinates of field of points.
    :return: list of list with the FOP after clustering.py.
    """
    db = DBSCAN(eps=0.5, min_samples=2).fit(fop)
    labels = db.labels_
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    #print(n_clusters_)
    longest = ""
    clusters = {}
    for i in range(n_clusters_):
        clusters["Cluster_"+str(i+1)] = []

    for i, coord in enumerate(fop):
        if labels[i] != -1:
            clusters["Cluster_"+str(labels[i]+1)].append(coord)
        else:
            pass

    for i in clusters:
        #print(len(clusters[i]))
        if longest == "" or len(clusters[i]) > len(clusters[longest]):
            longest = i
    #print(len(clusters[longest]))
    return clusters[longest]


def get_jaccard_distance(reference_fop, segment_fop, resolution):
    """
    Function that calculates the Jaccard distance between the points in reference_fop and the segment_fop. Uses the
    distance between points to calculate the intersection.

    :param reference_fop: list of list containing the reference FOP.
    :param segment_fop: list of list containing the segment FOP.
    :param resolution: resolution used to create FOP.
    :return: float with the Jaccard distance.
    """
    # Obtaining the trees for both field of points
    reference_tree = sp.spatial.cKDTree(reference_fop)
    segment_tree = sp.spatial.cKDTree(segment_fop)

    # Obtain the points that are at less than resolution/2.5 (aka have the same coordinates)
    clash_indices = reference_tree.query_ball_tree(segment_tree, resolution / 2.5, p=2, eps=0)

    # Count the points that intersect and convert to float
    intersection = len([x for x in clash_indices if x])/1.0

    # Obtain the union of both FOP
    union = float(len(reference_fop) + len(segment_fop) - intersection)

    # Calculate Jaccard distance
    jaccard = 1 - intersection / union

    return jaccard


def point_in_hull(point, hull, tolerance=1e-12):
    """
    a point is in the hull if and only if for every equation (describing the facets) the dot product between the point
    and the normal vector (eq[:-1]) plus the offset (eq[-1]) is less than or equal to zero. You may want to compare to a
    small, positive constant tolerance = 1e-12 rather than to zero because of issues of numerical precision (otherwise,
    you may find that a vertex of the convex hull is not in the convex hull).

    :param point:
    :param hull:
    :param tolerance:
    :return:
    """
    return all(
        (np.dot(eq[:-1], point) + eq[-1] <= tolerance)
        for eq in hull.equations)


def remove_convex_fop(trimmed_fop, trimmed_alpha):
    """

    :param trimmed_fop:
    :param trimmed_alpha:
    :return:
    """
    #todo add documentation
    points_in_hull = []
    trimmed_alpha_convex = sp.spatial.ConvexHull(trimmed_alpha)
    for point in trimmed_fop:
        if point_in_hull(point, trimmed_alpha_convex):
            points_in_hull.append(point)
    return points_in_hull


def calculate_pocket_gyration(pocket):
    """
    calculate_pocket_gyration is a function that takes the xyz coordinates of a field of points and calculates the
    radius of gyration. It assumes a mass of one for all the points.
    :param pocket: list of xyz coordinates of the field of points
    :return: float with the radius of gyration
    """
    def get_centroid(pocket):
        """
        Calculates the centroid or the center of mass for equal masses of a field of points.
        :param pocket: list of xyz coordinates of the field of points
        :return: list with coordinates od the centroid [x, y, z]
        """
        x = np.mean([x[0] for x in pocket])
        y = np.mean([y[1] for y in pocket])
        z = np.mean([z[2] for z in pocket])
        return [x, y, z]

    centroid = get_centroid(pocket)
    mass = len(pocket)
    gyration_radius = 0
    for i in pocket:
        gyration_radius += (calculate_distance_two_points(i, centroid))**2

    return np.sqrt(gyration_radius/mass)


if __name__ == "__main__":
    # Get arguments and load json settings file.
    parser = argparse.ArgumentParser(
        description="Obtain jaccard distance using a reference, a topology file and a MD trajectory file.")
    parser.add_argument("crd_file", type=str, help="Define the coordiante file. It is required")
    parser.add_argument("settings", type=str, help="Define the json file with the settings. It is required")
    parser.add_argument("--csv", type=str, help="will save results in an csv file with the provided filename")
    parser.add_argument("--we", action="store_true", help="Use this when you are using this for a progress coordinate "
                                                         "calculation in a weighted ensemble simulation")

    args = parser.parse_args()
    # Function that obtains and checks for settings in the settings json file.
    settings = check_input_settings_file(args.settings)
    # Checking the the trajectory doesn't have a rattle error.
    try:
        with open("seg.log", "r") as f:
            line = f.readline()
            while line:
                if "RATTLE" in line:
                    print("There was an error with the simulation")
                    sys.exit("There was an error with the simulation")
                else:
                    line = f.readline()
    except:
        print("Could not check log file. There could be a problem with the simulation")

    # Load reference pdb file and trajectory to then align the trajectory using the reference.
    ref_universe = MDAnalysis.Universe(settings["reference"])
    reference = ref_universe.select_atoms("protein")
    ensemble = MDAnalysis.Universe(settings["topology"], args.crd_file)

    # open file with selection string
    with open(settings["selection_file"], "r") as f:
        selection_pocket = f.readlines()[0]

    # open file with reference field of points
    if settings["fop_filetype"] == "xyz":
        reference_fop = parse_xyz_fop(settings["reference_fop"])
    elif settings["fop_filetype"] == "pdb":
        reference_fop = parse_pdb_fop(settings["reference_fop"]) #todo make parse_pdb_fop function
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

    # Define the number of times that the progress coordinates and auxiliary info will be calculated.
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
        if "bb_rmsd" not in results:
            results["bb_rmsd"] = []

    # this section is for the WESTPA analysis tools to work. The first point must be initial point.
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

        if "rog_pocket" in settings["auxdata"]:
            with open("parent_rog.txt", "r") as f:
                results["rog_pocket"].append(float(f.readlines()[-1]))

        if "bb_rmsd" in settings["auxdata"]:
            with open("parent_bb.txt", "r") as f:
                results["bb_rmsd"].append(float(f.readlines()[-1]))

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
            sigma = (1 - len(reference.select_atoms(selection_alignment))/len(reference.select_atoms("backbone"))) / 2

    # loop through frames of walker and calculate pcoord and auxdata
    for frame in np.linspace(0, len(ensemble.trajectory) - 1, num_points, dtype=int):
        ensemble.trajectory[frame]
        protein = ensemble.select_atoms("protein")
        # align the frame to the reference
        align.alignto(protein, reference, select=selection_alignment)
        # calculate pocket RMSD if needed
        if "prmsd" in results.keys():
            results["prmsd"].append(MDAnalysis.analysis.rms.rmsd(pocket_reference.positions,
                                                        ensemble.select_atoms(selection_pocket).positions))
        # The next lines calculate the Jaccard distance of the pocket comparing it to the reference
        if "jd" in results.keys() or "fops" in results.keys() or "pvol" in results.keys() or "rog_pocket" in results.keys():
            frame_coordinates = ensemble.select_atoms("protein").positions
            pocket_calpha = ensemble.select_atoms(selection_pocket + " and name CA*").positions
            frame_fop = get_field_of_points_dbscan(frame_coordinates, pocket_calpha, settings["center"],
                                        settings["resolution"], settings["radius"])
            results["jd"].append(get_jaccard_distance(reference_fop, frame_fop, settings["resolution"]))

        if "fops" in results.keys():
            results["fops"].append(frame_fop)

        if "pvol" in results.keys():
            results["pvol"].append(len(frame_fop) * (settings['resolution'] ** 3))

        if "rog_pocket" in results.keys():
            results["rog_pocket"].append(calculate_pocket_gyration(frame_fop))

        if "bb_rmsd" in results.keys():
            align.alignto(protein, reference, select="backbone")
            results["bb_rmsd"].append(MDAnalysis.analysis.rms.rmsd(reference.select_atoms("backbone").positions,
                                                                   protein.select_atoms("backbone").positions))

        if "composite" in results.keys():
            results["composite"].append(results["prmsd"][-1] + (sigma * results["bb_rmsd"][-1]))

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
            points_to_xyz_file("fop.txt", results["fops"], settings["resolution"], settings["radius"])
        elif settings["fop_filetype"] == "pdb":
            points_to_pdb("fop", results["fops"])
        # not sure pickles work
        elif settings["fop_filetype"] == "pickle":
            with open("fop.txt", "wb") as f:
                pickle.dump(frame_fop, f)

    if "pvol" in settings["auxdata"]:
        with open("pvol.txt", "w") as f:
            for i in results["pvol"]:
                f.write(str(i)+"\n")

    if "rog_pocket" in settings["auxdata"]:
        with open("rog.txt", "w") as f:
            for i in results["rog_pocket"]:
                f.write(str(i)+"\n")

    if "bb_rmsd" in settings["auxdata"]:
        with open("bb_rmsd.txt", "w") as f:
            for i in results["bb_rmsd"]:
                f.write(str(i)+"\n")

    if "prmsd" in settings["auxdata"]:
        with open("prmsd.txt", "w") as f:
            for i in results["prmsd"]:
                f.write(str(i)+"\n")

    if "jd" in settings["auxdata"]:
        with open("jd.txt", "w") as f:
            for i in results["jd"]:
                f.write(str(i)+"\n")

    if "composite" in settings["auxdata"]:
        with open("composite.txt", "w") as f:
            for i in results["composite"]:
                f.write(str(i)+"\n")

    if args.csv is not None:
        import pandas as pd
        results_pd = pd.DataFrame(results)
        results_pd.to_csv(args.csv)
