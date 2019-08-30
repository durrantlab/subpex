import argparse
import MDAnalysis
import json
import numpy as np
import scipy as sp
from MDAnalysis.analysis import rms, align


def get_rmsd_pocket(reference, protein, center, radius):
    """

    :param reference:
    :param protein:
    :param center:
    :param radius:
    :return:
    """
    # select heavy atoms within the sphere created by the center and the radius in the reference
    segment, transformation = prody.superpose(protein, reference)
    ref_pocket = reference.select("not water and not hydrogen and within " + str(radius) + " of pocketcenter",
                           pocketcenter=np.array(center))
    #print(type(ref_pocket))
    # Obtain the selection from above so we can pass it to the protein
    pocket_selection = ref_pocket.getSelstr()
    # Select heavy atoms in protein
    protein_pocket = segment.select(pocket_selection)
    #print(type(protein_pocket))
    # super pose reference and protein
    #protein_pocket = prody.superpose(protein_pocket, ref_pocket)
    # calculate rmsd
    pocket_rmsd = prody.calcRMSD(ref_pocket, protein_pocket)
    return pocket_rmsd


def points_to_pdb(filename, coordinates):
    """

    :param filename:
    :param coordinates:
    :return:
    """
    with open(filename, "w") as f:
        counter = 1
        for i in coordinates:
            text = "ATOM" + "{:>7}".format(counter)
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
            counter += 1
        f.write("\n")


def points_to_xyz_file(filename, coordinates, resolution, radius):
    """
    This Function takes the coordinates and the resolution and coordinates to write an xyz file that can be read by VMD.

    :param string filename: name for the xyz file to be created.
    :param list of lists coordinates: contains XYZ coordinates for each atom.
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

    :param list point1: coordinates of first point. [x, y, z]
    :param list point2: coordinates of second point. [x, y, z]
    :return float distance: distance between point1 and point2.
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
    VDW_radius = 1.4                                                                    #radius of water
    clash_indices = atoms_tree.query_ball_tree(fop_tree, VDW_radius, p=2, eps=0)
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
    # Count the points that intersect
    intersection = 0
    for index in clash_indices:
        for item in index:
            intersection += 1.0

    # Obtain the union of both FOP
    union = float(len(reference_fop) + len(segment_fop) - intersection)

    # Calculate Jaccard distance
    jaccard = 1 - intersection / union

    return jaccard


def get_field_of_points(protein, alphas, center, resolution, radius):
    """
    This function will take the coordinates of the protein and the coordinates of alpha carbons calculate the field of
    points (FOP).

    :param list of list protein: coordinates of protein atoms.
    :param list of lists alphas: coordinates of all alpha carbons
    :param center: center of the pocket (user defined).
    :param resolution: resolution for the FOP (user defined).
    :param radius: radius of sphere for FOP (user defined).
    :return list of lists: pocket shape as a FOP.
    """
    #get starting FOP and snapped center
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
    # a point is in the hull if and only if for every equation (describing the facets) the dot product between the point and
    # the normal vector (eq[:-1]) plus the offset (eq[-1]) is less than or equal to zero. You may want to compare to a small,
    # positive constant tolerance = 1e-12 rather than to zero because of issues of numerical precision (otherwise, you may
    # find that a vertex of the convex hull is not in the convex hull).
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
        if(point_in_hull(point, trimmed_alpha_convex)):
            points_in_hull.append(point)
    return points_in_hull


if __name__ == "__main__":
    # Get arguments and load json settings file.
    parser = argparse.ArgumentParser(description="Obtain jaccard distance using a reference, a topology file and a MD trajectory file.")
    parser.add_argument("reference", type=str, help="Define the reference PDB file. It is required")
    parser.add_argument("topology", type=str, help="Define the topology file (e.g. prmtop). It is required")
    parser.add_argument("dcd_file", type=str, help="Define the dcd file. It is required")
    parser.add_argument("settings", type=str, help="Define the json file with the settings. It is required")
    parser.add_argument("--csv", action="store_true", help="will save results in an csv file")
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
    ensemble = MDAnalysis.Universe(args.topology, args.dcd_file)
    alignment = align.AlignTraj(ensemble, reference, in_memory=True, select="backbone")
    alignment.run()

    # Define the pocket atoms as all atoms that are at a radius distance of the center point as defined by the user. 
    pocket_reference = reference.select_atoms("point {} {} {} {}".format(str(settings["center"][0]), 
        str(settings["center"][1]), str(settings["center"][2]), str(settings["radius"])))

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

    # Obtain RMSD of  the backbone PROBABLY WILL HAVE TO CHANGE to backbone atoms of the pocket.
    if settings["calculated_points"] == -1:
        num_points = len(ensemble.trajectory)
    else:
        num_points = settings["calculated_points"]

    jaccard = []
    pocket_rmsd = []
    bb_rmsd = []
    pvol = []

    for frame in np.linspace(0, len(ensemble.trajectory)-1, num_points, dtype=int):
        ensemble.trajectory[frame]
        bb_rmsd.append(MDAnalysis.analysis.rms.rmsd(reference.select_atoms("backbone").positions, ensemble.select_atoms("backbone").positions))
        pocket_rmsd.append(MDAnalysis.analysis.rms.rmsd(pocket_reference.positions, ensemble.select_atoms(selection_pocket).positions))
        frame_coordinates = ensemble.select_atoms("protein").positions
        frame_calpha = ensemble.select_atoms("name CA").positions
        frame_fop = get_field_of_points(frame_coordinates, frame_calpha, settings["center"],
                                          settings["resolution"], settings["radius"])

        pvol.append(len(frame_fop) * (settings['resolution'] ** 3))
        jaccard.append(get_jaccard_distance(reference_fop, frame_fop, settings["resolution"]))

    
    for i, jd in enumerate(jaccard):
        print(str(jd)+"    "+str(bb_rmsd[i]))


    if args.csv:
        import pandas as pd
        res_dict = {}
        res_dict["Jaccard"] = jaccard
        res_dict["Pocket RMSD"] = pocket_rmsd
        res_dict["Backbone RMSD"] = bb_rmsd
        res_dict["Pocket volume"] = pvol
        results = pd.DataFrame(res_dict)
        results.to_csv("results.csv")


