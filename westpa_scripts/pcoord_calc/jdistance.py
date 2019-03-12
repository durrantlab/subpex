import argparse
import prody
import json
import numpy as np
import scipy as sp

def points_to_xyz_file(filename, coordinates, resolution, dimension):
    volume = len(coordinates)*(resolution**3)
    with open(filename, 'w') as f:
        f.write(str(len(coordinates))+'\n')
        f.write('Generated in SubPEx. Res: {res}    Dim: {dim} Ang    Volume: {vol} Ang^3 \n'.format(res=resolution,
                                                                                        dim=dimension, vol = volume))
        for i in coordinates:
            f.write('CA    {}    {}    {} \n'.format(i[0], i[1], i[2]))


def calculate_distance_two_points(point1, point2):
    '''
    This function calculates the distance between two points in a 3D coordinate system.
    :param point1: coordinates of first point contained in a list. [x, y, z]
    :param point2: coordinates of second point contained in a list. [x, y, z]
    :return: distance between point1 and point2.
    '''
    distance = 0
    for i in range(3):
        distance += (point1[i]-point2[i])**2
        #print(np.sqrt(distance))
    return np.sqrt(distance)

def field_of_points(center, resolution, dimension):
    '''
    Function thathblahhhhh
    :return:
    '''
    # This part of the code will get a new center, a snapped center
    snapped_x = np.around(center[0] / resolution) * resolution
    snapped_y = np.around(center[1] / resolution) * resolution
    snapped_z = np.around(center[2] / resolution) * resolution
    center = np.around(np.array([snapped_x, snapped_y, snapped_z], float), 3)

    # Getting the field of points around the center.
    xs = np.round(np.arange(center[0] - dimension/2, center[0] + dimension/2 + resolution, resolution),3)
    ys = np.round(np.arange(center[1] - dimension/2, center[1] + dimension/2 + resolution, resolution),3)
    zs = np.round(np.arange(center[2] - dimension/2, center[2] + dimension/2 + resolution, resolution),3)
    field_of_points = np.vstack(np.meshgrid(xs,ys,zs)).reshape(3,-1).T
    # making the field of points a sphere
    field_of_points = [x for x in field_of_points if calculate_distance_two_points(x, center) < dimension / 2.0]

    return field_of_points, center


def get_trimmed_fop(field_of_points, atoms_points):
    #hgenerating cKDTrees to obtain distances
    atoms_tree = sp.spatial.cKDTree(atoms_points)
    fop_tree = sp.spatial.cKDTree(field_of_points)
    # get indices of fop_tree that are close to atoms tree
    clash_indices = atoms_tree.query_ball_tree(fop_tree, 1.5, p=2, eps=0) #TODO
    for index in clash_indices:
        for item in index:
            field_of_points[item] = None

    trimmed_fop = []
    for i in field_of_points:
        if type(i) == np.ndarray:
            trimmed_fop.append(i)
        else:
            pass

    return trimmed_fop


def get_jaccard_distance(reference_fop, pocket_fop):
    intersection = []
    union = reference_fop.copy()
    for i in pocket_fop:
        if i in union:
            intersection.append(i)
        else:
            union.append(i)

    return 1 - float(len(intersection)) / len(union)


if __name__ == "__main__":
    # Read the input parameter for the script given by the user.
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

    # Get the parameters to be used for calculating the Jaccard distance.
    protein_reference = prody.parsePDB(args.reference)
    protein = protein_reference.select('protein')
    ensemble = prody.parseDCD(args.dcd_file)
    n_frames = ensemble.numCoordsets()
    ensemble.setAtoms(protein_reference)
    ensemble.setCoords(protein_reference)
    ensemble.setAtoms(protein)
    ensemble.superpose()
    coords_protein =  ensemble.getCoords(n_frames)
    trimmed_coords_protein = [x for x in coords_protein if calculate_distance_two_points(x, settings['center']) < settings['dimension'] / 2.0]
    points_to_xyz_file('coords_prot.xyz', trimmed_coords_protein, settings['resolution'], settings['dimension'])

    fop, center = field_of_points(settings['center'], settings['resolution'], settings['dimension'])
    #points_to_xyz_file('fop.xyz', fop, settings['resolution'], settings['dimension'])
    trimmed_fop = get_trimmed_fop(fop, trimmed_coords_protein)
    points_to_xyz_file('pocket.xyz', trimmed_fop, settings['resolution'], settings['dimension'])

    ensemble.setAtoms(protein.calpha)
    coords_alpha = ensemble.getCoords(n_frames)
    trimmed_alpha = [x for x in coords_alpha if calculate_distance_two_points(x, settings['center']) < (settings['dimension']*1.25) / 2.0] #MAGIC
    points_to_xyz_file('alpha.xyz', trimmed_alpha, settings['resolution'], settings['dimension'])

