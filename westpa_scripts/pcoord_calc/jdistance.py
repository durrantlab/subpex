import argparse
import prody
import json
import numpy as np
import scipy as sp


def points_to_xyz_file(filename, coordinates, resolution, radius):
    '''
    This Function takes the coordinates and the resolution and coordinates to write an xyz file that can be read by VMD.

    :param string filename: name for the xyz file to be created.
    :param list of lists coordinates: contains XYZ coordinates for each atom.
    :param float resolution: Resolution in Angstroms.
    :param float radius: radius in Angstrom for FOP.
    '''
    # calculate the volume using the resolution
    volume = len(coordinates) * (resolution ** 3)

    # Writing the xyz file.
    with open(filename, 'w') as f:
        f.write(str(len(coordinates)) + '\n')
        f.write('Generated in SubPEx. Res: {res}    Rad: {dim} Ang    Volume: {vol} Ang^3 \n'.format(res=resolution,
                                                                                                     dim=radius,
                                                                                                     vol=volume))
        # Adding each coordinate as a C-alpha for visualization purposes
        for i in coordinates:
            f.write('CA    {}    {}    {} \n'.format(i[0], i[1], i[2]))


def calculate_distance_two_points(point1, point2):
    '''
    This function calculates the distance between two points in a 3D coordinate system.

    :param list point1: coordinates of first point. [x, y, z]
    :param list point2: coordinates of second point. [x, y, z]
    :return float distance: distance between point1 and point2.
    '''
    distance = 0
    for i in range(3):
        distance += (point1[i] - point2[i]) ** 2

    return np.sqrt(distance)


def field_of_points(center, resolution, radius):
    '''
    Function that snaps the provided center to the a center congruent to the resolution and radiuss given by the
    user. Then it generates a spherical field of points (FOP) with radius = radius/2 with points at with a distance
    between them equal to the resolution.

    :param list center: provided center coordinates. [X, Y, Z]
    :param float resolution: Distance between points which give resolution of the pocket volume calculator method.
    :param float radius: Defines the radius of the sphere of the FOP.
    :return list of lists : a spherical FOP
    '''
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
    '''
    Removes points that clash with atoms within the FOP.

    :param field_of_points: List of list containing all the points for hte FOP.
    :param atoms_points: List of list containing the coordinates of atoms within the FOP.
    :return: Returns a list of list that has the trimmed FOP.
    '''
    # Generating cKDTrees to obtain distances
    atoms_tree = sp.spatial.cKDTree(atoms_points)
    fop_tree = sp.spatial.cKDTree(field_of_points)

    # get indices of fop_tree that are close to atoms_tree
    VDW_radius = 1.4                                                                    #radius of water
    clash_indices = atoms_tree.query_ball_tree(fop_tree, VDW_radius, p=2, eps=0)
    for index in clash_indices:
        for item in index:
            field_of_points[item] = None

    # Creating the FOP that does not contain points that clash with the atoms provided tpo the function.
    trimmed_fop = []
    for i in field_of_points:
        if type(i) == np.ndarray:
            trimmed_fop.append(i)
        else:
            pass

    return trimmed_fop


def get_jaccard_distance(reference_fop, pocket_fop):
    '''
    Function that calculates the Jaccard distance between the points in reference_fop and the pocket_fop.

    :param reference_fop: list of list containing the reference FOP.
    :param pocket_fop: list of list containing the pocket FOP.
    :return: float with the Jaccard distance.
    '''
    intersection = []
    union = reference_fop.copy()
    for i in pocket_fop:
        if i in union:
            intersection.append(i)
        else:
            union.append(i)

    return 1 - float(len(intersection)) / len(union)


def get_field_of_points(protein, alphas, center, resolution, radius):
    '''
    This function will take the coordinates of the protein and the coordinates of alpha carbons calculate the field of
    points (FOP).

    :param list of list protein: coordinates of protein atoms.
    :param list of lists alphas: coordinates of all alpha carbons
    :param center: center of the pocket (user defined).
    :param resolution: resolution for the FOP (user defined).
    :param radius: radius of sphere for FOP (user defined).
    :return list of lists: pocket shape as a FOP.
    '''
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


def remove_convex_fop(trimmed_fop, trimmed_alpha):
    '''

    :param trimmed_fop:
    :param trimmed_alpha:
    :return:
    '''
    return trimmed_fop


if __name__ == "__main__":
    # Get arguments and load json settings file.
    parser = argparse.ArgumentParser(description="Obtain jaccard distance using a reference and a MD trajectory file.")
    parser.add_argument("reference", type=str, help="Define the reference PDB file. It is required")
    parser.add_argument("dcd_file", type=str, help="Define the dcd file. It is required")
    parser.add_argument("settings", type=str, help="Define the json file with the settings. It is required")
    parser.add_argument("--debug", action="store_true", help="To keep all the files to debug any problem")
    args = parser.parse_args()

    try:
        with open(args.settings, "r") as f:
            settings = json.load(f)
    except IOError:
        print("Could not load the json file with the settings")
        print("make sure the file exists and is correctly formatted")
        raise IOError("Could not load the json file with the setttings")

    # Load pdb file which should be the reference and obtained the pocket field of points.
    protein_reference = prody.parsePDB(args.reference)
    protein = protein_reference.select('protein')
    reference_coordinates = protein.getCoords()
    reference_alpha = protein.calpha.getCoords()
    reference_fop = get_field_of_points(reference_coordinates, reference_alpha, settings['center'],
                                        settings['resolution'], settings['radius'])
    points_to_xyz_file('ref_pocket.xyz', reference_fop, settings['resolution'], settings['radius'])

    # Load dcd file and grab last frame to obtained the pocket field of points of that frames conformation.
    ensemble = prody.parseDCD(args.dcd_file)
    n_frames = ensemble.numCoordsets()
    ensemble.setAtoms(protein_reference)
    ensemble.setCoords(protein_reference)
    ensemble.setAtoms(protein)
    ensemble.superpose()
    coords_protein = ensemble.getCoords(n_frames)
    ensemble.setAtoms(protein.calpha)
    coords_alpha = ensemble.getCoords(n_frames)
    segment_fop = get_field_of_points(reference_coordinates, coords_alpha, settings['center'], settings['resolution'],
                                       settings['radius'])

    points_to_xyz_file('seg_pocket.xyz', segment_fop, settings['resolution'], settings['radius'])
    #points_to_xyz_file('alpha.xyz', trimmed_alpha, settings['resolution'], settings['radius'])

    # Compare pockets of reference and the segment and calculate Jaccard distance between both.
    #jaccard = get_jaccard_distance(reference_fop, segment_fop)
