import argparse
import prody
import json
import numpy as np
import scipy as sp

def get_rmsd_pocket(reference, protein, center, radius):
    '''

    :param reference:
    :param protein:
    :param center:
    :param radius:
    :return:
    '''

    # Align the protein pocket onto the reference pocket.
    protein_pocket, ref_pocket = align_last_frame_to_ref_by_pocket(protein, reference, center, radius, False)

    # calculate rmsd
    pocket_rmsd = prody.calcRMSD(ref_pocket, protein_pocket)

    return pocket_rmsd


def points_to_pdb(filename, coordinates):
    '''

    :param filename:
    :param coordinates:
    :return:
    '''
    with open(filename, 'w') as f:
        counter = 1
        for i in coordinates:
            text = 'ATOM' + '{:>7}'.format(counter)
            text += '{:^6}'.format('CA')
            text += '{:>3}'.format('ALA')
            text += '{:>6}'.format('A')
            text += '{:>12}'.format(i[0])
            text += '{:>8}'.format(i[1])
            text += '{:>8}'.format(i[2])
            text += '{:>6}'.format('1.0')
            text += '{:>6}'.format('1.0')
            text += '\n'
            f.write(text)
            counter += 1
        f.write('\n')


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
    if type(atoms_points) is list and len(atoms_points) > 0:
        atoms_points = np.vstack(atoms_points)
    if type(field_of_points) is list and len(field_of_points) > 0:
        field_of_points = np.vstack(field_of_points)

    #if len(atoms_points) > 0 or len(field_of_points) > 0:
    #    return 1.0

    atoms_tree = sp.spatial.cKDTree(atoms_points)
    fop_tree = sp.spatial.cKDTree(field_of_points)

    # get indices of fop_tree that are close to atoms_tree
    #VDW_radius = 1.4                                                                    #radius of water
    VDW_radius = 2.0

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


def get_jaccard_distance(reference_fop, last_frame_fop, resolution):
    '''
    Function that calculates the Jaccard distance between the points in reference_fop and the last_frame_fop. Uses the
    distance between points to calculate the intersection.

    :param reference_fop: list of list containing the reference FOP.
    :param last_frame_fop: list of list containing the segment FOP.
    :param resolution: resolution used to create FOP.
    :return: float with the Jaccard distance.
    '''
    # Obtaining the trees for both field of points
    reference_tree = sp.spatial.cKDTree(reference_fop)
    segment_tree = sp.spatial.cKDTree(last_frame_fop)
    # Obtain the points that are at less than resolution/2.5 (aka have the same coordinates)
    clash_indices = reference_tree.query_ball_tree(segment_tree, resolution / 2.5, p=2, eps=0)
    # Count the points that intersect
    intersection = 0
    for index in clash_indices:
        for item in index:
            intersection += 1.0

    # Obtain the union of both FOP
    union = float(len(reference_fop) + len(last_frame_fop) - intersection)

    # Calculate Jaccard distance
    jaccard = 1 - intersection / union

    return jaccard


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


def point_in_hull(point, hull, tolerance=1e-12):
    # a point is in the hull if and only if for every equation (describing the facets) the dot product between the point and
    # the normal vector (eq[:-1]) plus the offset (eq[-1]) is less than or equal to zero. You may want to compare to a small,
    # positive constant tolerance = 1e-12 rather than to zero because of issues of numerical precision (otherwise, you may
    # find that a vertex of the convex hull is not in the convex hull).
    return all(
        (np.dot(eq[:-1], point) + eq[-1] <= tolerance)
        for eq in hull.equations)


def remove_convex_fop(trimmed_fop, trimmed_alpha):
    '''

    :param trimmed_fop:
    :param trimmed_alpha:
    :return:
    '''
    points_in_hull = []
    trimmed_alpha_convex = sp.spatial.ConvexHull(trimmed_alpha)
    for point in trimmed_fop:
        if(point_in_hull(point, trimmed_alpha_convex)):
            points_in_hull.append(point)
    return points_in_hull

def align_last_frame_to_ref_by_pocket(last_frame, ref_protein, center, radius, return_whole_protein):
    # select heavy atoms within the sphere created by the center and the radius in the reference
    ref_pocket = ref_protein.select('not water and not hydrogen and not resname "Cl-" "Na+" SOD CLA and within ' + str(radius) + ' of pocketcenter',
                                    pocketcenter=np.array(center))

    # Obtain the selection from above so we can pass it to the protein
    ref_pocket_sel = ref_pocket.getSelstr()

    # Select heavy atoms in protein
    last_frame_pocket = last_frame.select(ref_pocket_sel)

    # super pose the last-frame pocket onto the reference pocket. Interestingly, this also
    # aligns last_frame.
    last_frame_pocket, transform = prody.superpose(last_frame_pocket, ref_pocket)

    if return_whole_protein:
        # Apply that transform to the whole protein and return that.
        # Nevermind. last_frame is also aligned when you align by the
        # pocket above.
        # last_frame = transform.apply(last_frame)
        return last_frame, ref_protein
    else:
        # Return only the aligned pocket atoms from the selection. For RMSD calc.
        return last_frame_pocket, ref_pocket

if __name__ == "__main__":

    DEBUG = True

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
        raise IOError("Could not load the json file with the settings")

    # Load pdb file and dcd file to load the reference and the trajectory to obtain field of points of the pocket.
    protein = prody.parsePDB(args.reference)
    ref_protein = protein.select('protein')
    ensemble = prody.parseDCD(args.dcd_file)
    ensemble.setAtoms(ref_protein)

    # Save the last frame of the trajectory file as a PDB file and then load it. This is done because the way
    # ensembles work does not let us superpose and change the coordinates of the frame.
    prody.writePDB('seg.last_frame.aligned_to_ref_pocket.pdb', ensemble[-1])
    last_frame = prody.parsePDB('seg.last_frame.aligned_to_ref_pocket.pdb')

    # Superpose the last_frame of the trajecotry file and move the coordinate system of the last frame.
    # Below aligns by the whol protein. Commented out, because I want to align by pocket here too.
    #last_frame, _transform = prody.superpose(last_frame, ref_protein)
    last_frame, ref_protein = align_last_frame_to_ref_by_pocket(last_frame, ref_protein, settings['center'], settings['radius'], True)

    # Write the nwe protein in the new coordinate system.
    prody.writePDB('seg.last_frame.aligned_to_ref_pocket.pdb', last_frame)

    if DEBUG: prody.writePDB("ref_protein.pdb", ref_protein)

    # Obtain RMSD of the backbone PROBABLY WILL HAVE TO CHANGE to backbone atoms of the pocket.
    # rmsd_all_backbone never used.
    # rmsd_all_backbone = prody.calcRMSD(ref_protein.backbone, last_frame.backbone)

    # Obtain coordinates for reference atoms and coordinates for alpha carbons.
    # JDD comment: Really necssary to calculate this every time? Probably fast anyway.
    ref_coors = ref_protein.getCoords()
    ref_alpha = ref_protein.calpha.getCoords()
    reference_fop = get_field_of_points(ref_coors, ref_alpha, settings['center'],
                                        settings['resolution'], settings['radius'])

    # Save xyz coordinates of FOP
    if DEBUG: points_to_xyz_file('ref_pocket.xyz', reference_fop, settings['resolution'], settings['radius'])

    # Obtain coordinates for last frame atoms and coordinates for alpha carbons.
    last_frame_coors = last_frame.getCoords()
    last_frame_alpha = last_frame.calpha.getCoords()
    last_frame_fop = get_field_of_points(last_frame_coors, last_frame_alpha, settings['center'],
                                         settings['resolution'], settings['radius'])

    if DEBUG: points_to_xyz_file('seg_pocket.xyz', last_frame_fop, settings['resolution'], settings['radius'])

    # Compare pockets of reference and the segment and calculate Jaccard distance between both.
    jaccard = get_jaccard_distance(reference_fop, last_frame_fop, settings['resolution'])

    # Calcult the rmsd too.
    rmsd = get_rmsd_pocket(ref_protein, last_frame, settings['center'], settings['radius'])

    # Print both metrics to the screen.
    pcoord = str(jaccard)+'    '+str(rmsd)
    print(pcoord)

    # Save data (for debugging?)
    with open('jaccard.dat', 'w') as f:
        f.write('1 '+str(jaccard))

    with open('rmsd.dat', 'w') as f:
        f.write('1 '+str(rmsd))

    with open('pvol.txt', 'a') as f:
        f.write(str(len(last_frame_fop) * (settings['resolution'] ** 3))+'\n')

    with open('timer.txt', 'w') as f:
        f.write('Done with the script')

