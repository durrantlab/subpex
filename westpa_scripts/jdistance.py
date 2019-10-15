import argparse
import MDAnalysis
import json
import numpy as np
import scipy as sp
from MDAnalysis.analysis import rms, align
import sys


def points_to_pdb(filename, coordinates):
    """
    points_to_pdb will write a pdb file full of C-alphas to be able to load the file in a visualisation software and
    represent the field of points.

    :param filename: (str) name of the pdb file to create.
    :param coordinates: (list of lists) XYZ coordinates of the field of points to write into pdb file.
    """
    with open(filename, "w") as f:
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
    volume = len(coordinates) * (resolution ** 3)

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
    for i in range(3):
        distance += (point1[i] - point2[i]) ** 2

    return np.sqrt(distance)


def field_of_points(center, resolution, radius):
    """
    Function that snaps the provided center to the a center congruent to the resolution and radiuss given by the
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


def get_field_of_points(protein, alphas, center, resolution, radius):
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
    # trim protein alpha atoms to those in the pocket.
    trimmed_alpha = [x for x in alphas if calculate_distance_two_points(x, center) < (radius * 1.25)]
    # remove points outside the convex hull created by the alpha carbons in the pocket.
    pocket = remove_convex_fop(trimmed_fop, trimmed_alpha)

    return pocket


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
    parser.add_argument("reference", type=str, help="Define the reference PDB file. It is required")
    parser.add_argument("topology", type=str, help="Define the topology file (e.g. prmtop). It is required")
    parser.add_argument("dcd_file", type=str, help="Define the dcd file. It is required")
    parser.add_argument("settings", type=str, help="Define the json file with the settings. It is required")
    parser.add_argument("--csv", type=str, help="will save results in an csv file with the provided filename")
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
    ref_universe = MDAnalysis.Universe(args.reference)
    reference = ref_universe.select_atoms("protein")
    ensemble = MDAnalysis.Universe(args.topology, args.dcd_file)

    # Define the pocket atoms as all atoms that are at a radius distance of the center point as defined by the user.
    pocket_reference = reference.select_atoms("point {} {} {} {}".format(str(settings["center"][0]),
                                                                         str(settings["center"][1]),
                                                                         str(settings["center"][2]),
                                                                         str(settings["radius"])))

    # select heavy atoms
    pocket_reference = pocket_reference.select_atoms("not name H")

    # Get the indexes of atoms selected in the reference
    indexes = pocket_reference.ix

    # Use the indexes to create a selection string to pass selection from reference to ensemble
    selection_pocket = "bynum " + str(indexes[0] + 1)
    for i in indexes[1:]:
        selection_pocket += " or bynum " + str(i + 1)

    # Obtain coordinates for reference atoms and coordinates for alpha carbons.
    reference_coordinates = reference.positions
    reference_alpha = reference.select_atoms("name CA").positions
    reference_fop = get_field_of_points(reference_coordinates, reference_alpha, settings["center"],
                                        settings["resolution"], settings["radius"]) #todo add function to trim the non contiguous points in the pocket.

    # Define the number of times that the progress coordinates and auxiliary info will be calculated.
    if settings["calculated_points"] == -1:
        num_points = len(ensemble.trajectory)
    else:
        num_points = settings["calculated_points"]

    # Create a list of all elements to calculate
    results = {}
    results["jaccard"] = []
    results["pocket_rmsd"] = []
    if args.csv != None:
        results["bb_rmsd"] = []
        results["pvol"] = []
        results["rog_pocket"] = []
    else:
        if args.pvol:
            results["pvol"] = []
        if args.rog:
            results["rog_pocket"] = []
        if args.bb_rmsd:
            results["bb_rmsd"] = []

    # go the the frames and do the calculations
    for frame in np.linspace(0, len(ensemble.trajectory) - 1, num_points, dtype=int):
        ensemble.trajectory[frame]
        protein = ensemble.select_atoms("protein")
        # align the frame to the reference
        align.alignto(protein, reference, select="backbone")
        # calculate pocket RMSD
        results["pocket_rmsd"].append(MDAnalysis.analysis.rms.rmsd(pocket_reference.positions,
                                                        ensemble.select_atoms(selection_pocket).positions))
        # The next lines calculate the jaccard distance of the pocket
        frame_coordinates = ensemble.select_atoms("protein").positions
        frame_calpha = ensemble.select_atoms("name CA").positions
        frame_fop = get_field_of_points(frame_coordinates, frame_calpha, settings["center"],
                                        settings["resolution"], settings["radius"])
        results["jaccard"].append(get_jaccard_distance(reference_fop, frame_fop, settings["resolution"]))

        # calculate backbone RMSD
        if args.bb_rmsd or (args.csv is not None):
            results["bb_rmsd"].append(MDAnalysis.analysis.rms.rmsd(reference.select_atoms("backbone").positions,
                                                    protein.select_atoms("backbone").positions))
        else:
            pass
        if args.pvol or (args.csv is not None):
            results["pvol"].append(len(frame_fop) * (settings['resolution'] ** 3))
        else:
            pass
        # calculate the radius of gyration of the pocket
        if args.rog or (args.csv is not None):
            results["rog_pocket"].append(calculate_pocket_gyration(frame_fop))
        else:
            pass

    # saving all the files required for westpa to run
    with open("pcoord.txt", "w") as f:
        for i, value in enumerate(results["jaccard"]):
            f.write("{:.4f}    {:.4f}\n".format(value, results["pocket_rmsd"][i]))

    if args.pvol:
        with open("pvol.txt", "w") as f:
            for i in results["pvol"]:
                f.write(str(i)+"\n")

    if args.rog:
        with open("rog.txt", "w") as f:
            for i in results["rog_pocket"]:
                f.write(str(i)+"\n")

    if args.bb_rmsd:
        with open("bb_rmsd.txt", "w") as f:
            for i in results["bb_rmsd"]:
                f.write(str(i)+"\n")

    if args.csv is not None:
        import pandas as pd
        results_pd = pd.DataFrame(results)
        results_pd.to_csv(args.csv)
