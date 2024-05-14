from typing import Any, List, Tuple

import argparse
import sys

import MDAnalysis
import numpy as np
import scipy as sp
from MDAnalysis.analysis import align
from sklearn.cluster import DBSCAN


def points_to_pdb(file_path: str, coordinates: List[List[float]]) -> None:
    """Writes a pdb file full of C-alphas so users can visualize the field of
    points in software such as VMD.

    Args:
        filename: name of the pdb file to create.
        coordinates: XYZ coordinates of the field of points to write to pdb file.
    """
    # TODO: Erich - modify to be able to do multiframe pdb
    with open(file_path, "w", encoding="utf-8") as f:
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


def parse_pdb_fop(filename: str) -> List[List[float]]:
    """Opens the field of points that is stored in a pdb file with name
    'filename', parse it and return the fop as a list of lists.

    Args:
        filename: name of the pdb file with the field of points.

    Returns:
        fop: XYZ coordinates of the field of points.
    """

    with open(filename, "r", encoding="utf-8") as f:
        text = f.readlines()
    fop = []
    for line in text:
        if line.startswith("ATOM"):
            split = line.split()
            x = float(split[4])
            y = float(split[5])
            z = float(split[6])
            fop.append([x, y, z])
        else:
            pass
    return fop


def parse_xyz_fop(filename: str) -> List[List[float]]:
    """Reads a filename and parses the field of points to convert it into a
    list of lists.

    Args:
        filename (str): filename of the field of points to open and parse

    Returns:
        fop (list of lists): x, y, z coordinates of the field of points
    """

    # open reference fop xyz file
    with open(filename, "r", encoding="utf-8") as f:
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


def points_to_xyz_file(
    filename: str, coordinates: List[List[float]], resolution: float, radius: float
) -> None:
    """Takes the coordinates and the resolution, writes write an xyz file that
    users can load into visualization software such as VMD.

    Args:
        filename: name for the xyz file to be created.
        coordinates: contains XYZ coordinates for each atom.
        resolution: resolution in Angstroms.
        radius: radius in Angstrom for the field of points.
    """

    # calculate the volume using the resolution
    # TODO: Erich - modify to be able to do multiframe xyz
    volume = len(coordinates) * (resolution**3)

    # Writing the xyz file.
    with open(filename, "w", encoding="utf-8") as f:
        f.write(str(len(coordinates)) + "\n")
        f.write(
            "Generated in SubPEx. Res: {res}    Rad: {dim} Ang    Volume: {vol} Ang^3 \n".format(
                res=resolution, dim=radius, vol=volume
            )
        )
        # Adding each coordinate as a C-alpha for visualization purposes
        for i in coordinates:
            f.write(f"CA    {i[0]}    {i[1]}    {i[2]} \n")


def calculate_distance_two_points(point1: List[float], point2: List[float]) -> float:
    """Calculates the distance between two points in a 3D coordinate system.

    Args:
        point1 (list): coordinates of first point. [x, y, z]
        point2 (list): coordinates of second point. [x, y, z]

    Returns:
        distance (float): distance between point1 and point2.
    """

    distance = 0.0
    for p1, p2 in zip(point1, point2):
        distance += (p1 - p2) ** 2

    return float(np.sqrt(distance))


def field_of_points(
    center: List[float], resolution: float, radius: float
) -> Tuple[List[List[float]], List[float]]:
    """Snaps the provided center to the a center congruent to the resolution
    and radius given by the user. Then generates a spherical field of points
    (FOP) with radius = radius/2 with points spaced "resolution" apart.

    Args:
        center (list): provided center coordinates. [X, Y, Z]
        resolution (float): Distance between points which give resolution of the
            pocket volume calculator method.
        radius (float): Defines the radius of the sphere of the FOP.

    Returns:
        field_of_points (list of lists): a spherical field of points
        center (list): center of the field of points
    """

    # This part of the code will get a new center, a snapped center.
    snapped_x = np.around(center[0] / resolution) * resolution
    snapped_y = np.around(center[1] / resolution) * resolution
    snapped_z = np.around(center[2] / resolution) * resolution
    center = np.around(np.array([snapped_x, snapped_y, snapped_z], float), 3)

    # Getting the field of points around the center.
    xs = np.round(
        np.arange(center[0] - radius, center[0] + radius + resolution, resolution), 3
    )
    ys = np.round(
        np.arange(center[1] - radius, center[1] + radius + resolution, resolution), 3
    )
    zs = np.round(
        np.arange(center[2] - radius, center[2] + radius + resolution, resolution), 3
    )
    field_of_points = np.vstack(np.meshgrid(xs, ys, zs)).reshape(3, -1).T

    # making the field of points a sphere of radius = radius/2.
    field_of_points = [
        x for x in field_of_points if calculate_distance_two_points(x, center) < radius
    ]

    return field_of_points, center


def get_trimmed_fop(
    field_of_points: List[List[float]], atoms_points: List[List[float]]
) -> List[List[float]]:
    """Removes points that clash with atoms within the field of points.

    Args:
        field_of_points (list of lists): a field of points to be trimmed of
            collisions with atoms
        atoms_points (list of lists): coordinates of atoms within the field of points.

    Returns:
        trimmed_fop (list of lists): a field of points without clashes
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

    # Creating the FOP that does not contain points that clash with the atoms
    # provided to the function.
    trimmed_fop = []
    for i in field_of_points:
        if type(i) == np.ndarray:
            trimmed_fop.append(i)
        else:
            pass

    return trimmed_fop


def get_field_of_points_dbscan(
    protein: List[List[float]],
    alphas: List[List[float]],
    center: List[float],
    resolution: float,
    radius: float,
) -> List[List[float]]:
    """Takes the coordinates of the protein and the coordinates of alpha
    carbons and calculates the field of points (FOP).

    Args:
        protein (list of lists): coordinates of protein atoms.
        alphas (list of lists): coordinates of all alpha carbons.
        center (list): center of the pocket (user defined).
        resolution (float): resolution for the FOP (user defined).
        radius (float): radius of sphere for FOP (user defined).

    Returns:
        pocket shape as a field of points.
    """

    # get starting FOP and snapped center
    fop, center = field_of_points(center, resolution, radius)
    # trim protein atoms to those within the FOP sphere.
    trimmed_coords_protein = [
        x for x in protein if calculate_distance_two_points(x, center) < radius
    ]
    # remove points in FOP that have steric clashes with protein.
    trimmed_fop = get_trimmed_fop(fop, trimmed_coords_protein)
    # remove points outside the convex hull created by the alpha carbons in the pocket.
    pocket = remove_convex_fop(trimmed_fop, alphas)
    pocket = cluster_dbscan(pocket)

    return pocket


def cluster_dbscan(fop: List[List[float]]) -> List[List[float]]:
    """Takes the coordinates of the field of points (fop) and uses a clustering
    algorithm to trim the fop to the points. Keeps the cluster that contains
    the center of the pocket.

    Args:
        fop: coordinates of field of points.

    Returns:
        field of points after clustering
    """

    if len(fop) == 0:
        # sometimes no points are present in the FOP.
        return []

    db = DBSCAN(eps=0.5, min_samples=2).fit(fop)
    labels = db.labels_
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

    if n_clusters_ == 0:
        # In some cases, even when fop has a few points, there are no clusters.
        # In this case, just return an empty list.
        return []

    longest = ""
    clusters: dict[str, List[Any]] = {}
    for i in range(n_clusters_):
        clusters["Cluster_" + str(i + 1)] = []

    for i, coord in enumerate(fop):
        if labels[i] != -1:
            clusters["Cluster_" + str(labels[i] + 1)].append(coord)
        else:
            pass

    for i in clusters:
        if longest == "" or len(clusters[i]) > len(clusters[longest]):
            longest = i

    return clusters[longest]


def get_jaccard_distance(
    reference_fop: List[List[float]], segment_fop: List[List[float]], resolution: float
) -> float:
    """Calculates the Jaccard distance between the points in reference_fop and
    the segment_fop. Uses the distance between points to calculate the intersection.

    Args:
        reference_fop: reference FOP.
        segment_fop: segment FOP.
        resolution: resolution used to create FOP.

    Returns:
        Jaccard distance.
    """

    if len(segment_fop) == 0:
        # sometimes no points are present in the FOP.
        return 1.0

    # Obtaining the trees for both field of points
    reference_tree = sp.spatial.cKDTree(reference_fop)
    segment_tree = sp.spatial.cKDTree(segment_fop)

    # Obtain the points that are at less than resolution/2.5 (aka have the same
    # coordinates)
    clash_indices = reference_tree.query_ball_tree(
        segment_tree, resolution / 2.5, p=2, eps=0
    )

    # Count the points that intersect and convert to float
    intersection = len([x for x in clash_indices if x]) / 1.0

    # Obtain the union of both FOP
    union = float(len(reference_fop) + len(segment_fop) - intersection)

    # Calculate Jaccard distance
    jaccard = 1 - intersection / union

    return jaccard


def point_in_hull(
    point: List[float], hull: scipy.spatial.ConvexHull, tolerance: float = 1e-12
) -> bool:
    """Determines if a point is contained in a hull. A point is in the hull if
    and only if for every equation (describing the facets), the dot product
    between the point and the normal vector (eq[:-1]) plus the offset (eq[-1])
    is less than or equal to zero. You may want to compare to a small, positive
    constant tolerance = 1e-12 rather than to zero because of issues of
    numerical precision (otherwise, you may find that a vertex of the convex
    hull is not in the convex hull).

    Args:
        point: coordinates of the point to be considered.
        hull: Convex hulls in N dimensions.
        tolerance: tolerance for the calculation.

    Returns:
        Returns True if point in hull, False otherwise
    """
    return all((np.dot(eq[:-1], point) + eq[-1] <= tolerance) for eq in hull.equations)


def remove_convex_fop(
    trimmed_fop: List[List[float]], trimmed_alpha: List[List[float]]
) -> List[List[float]]:
    """Uses a convex hull to trim points from the field of points that are
    outside of the protein pocket.

    Args:
        trimmed_fop: field of points trimmed of clashes with atoms.
        trimmed_alpha: alpha atoms to define the convex hull.

    Returns:
        points_in_hull: FOP trimmed of points outside the convex hull defined by the
            protein.
    """

    points_in_hull = []
    trimmed_alpha_convex = sp.spatial.ConvexHull(trimmed_alpha)
    for point in trimmed_fop:
        if point_in_hull(point, trimmed_alpha_convex):
            points_in_hull.append(point)
    return points_in_hull


def calculate_pocket_gyration(pocket: List[List[float]]) -> float:
    """Takes the xyz coordinates of a field of points and calculates the radius
    of gyration. It assumes a mass of one for all the points.

    Args:
        pocket: the field of points defining the pocket shape.

    Returns:
        radius of gyration of the pocket.
    """

    mass = len(pocket)

    if mass == 0:
        # Sometimes there are no points.
        return -9999

    centroid = get_centroid(pocket)
    gyration_radius = 0.0
    for i in pocket:
        gyration_radius += (calculate_distance_two_points(i, centroid)) ** 2

    return float(np.sqrt(gyration_radius / mass))


def get_centroid(pocket: List[List[float]]) -> List[float]:
    """Calculates the centroid or the center of mass for equal masses of a
    field of points.

    Args:
        pocket (list of lists):  the field of points defining the pocket shape.

    Returns:
        centroid (list): coordinates of the centroid of the pocket.
    """

    x = float(np.mean([x[0] for x in pocket]))
    y = float(np.mean([y[1] for y in pocket]))
    z = float(np.mean([z[2] for z in pocket]))
    return [x, y, z]


if __name__ == "__main__":
    # Get arguments and load json settings file.
    parser = argparse.ArgumentParser(
        description="Obtain jaccard distance using a reference, a topology file and a MD trajectory file."
    )
    parser.add_argument(
        "crd_file", type=str, help="Define the coordiante file. It is required"
    )
    parser.add_argument(
        "settings",
        type=str,
        help="Define the json file with the settings. It is required",
    )
    parser.add_argument(
        "--csv",
        type=str,
        help="will save results in an csv file with the provided filename",
    )
    parser.add_argument(
        "--we",
        action="store_true",
        help="Use this when you are using this for a progress coordinate "
        "calculation in a weighted ensemble simulation",
    )

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

    # Load reference pdb file and trajectory to then align the trajectory using
    # the reference.
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
        # TODO: Erich - make parse_pdb_fop function
        reference_fop = parse_pdb_fop(settings["reference_fop"])
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

    # Define the number of times that the progress coordinates and auxiliary
    # info will be calculated.
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
        if "bb" not in results:
            results["bb"] = []

    # this section is for the WESTPA analysis tools to work. The first point
    # must be initial point.
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

        if "rog" in settings["auxdata"]:
            with open("parent_rog.txt", "r") as f:
                results["rog"].append(float(f.readlines()[-1]))

        if "bb" in settings["auxdata"]:
            with open("parent_bb.txt", "r") as f:
                results["bb"].append(float(f.readlines()[-1]))

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
            sigma = (
                1
                - len(reference.select_atoms(selection_alignment))
                / len(reference.select_atoms("backbone"))
            ) / 2

    # loop through frames of walker and calculate pcoord and auxdata
    for frame in np.linspace(0, len(ensemble.trajectory) - 1, num_points, dtype=int):
        ensemble.trajectory[frame]
        protein = ensemble.select_atoms("protein")
        # align the frame to the reference
        align.alignto(protein, reference, select=selection_alignment)
        # calculate pocket RMSD if needed
        if "prmsd" in results.keys():
            results["prmsd"].append(
                MDAnalysis.analysis.rms.rmsd(
                    pocket_reference.positions,
                    ensemble.select_atoms(selection_pocket).positions,
                )
            )

        # The next lines calculate the Jaccard distance of the pocket comparing
        # it to the reference
        if (
            "jd" in results.keys()
            or "fops" in results.keys()
            or "pvol" in results.keys()
            or "rog" in results.keys()
        ):
            frame_coordinates = ensemble.select_atoms("protein").positions
            pocket_calpha = ensemble.select_atoms(
                selection_pocket + " and name CA*"
            ).positions
            frame_fop = get_field_of_points_dbscan(
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
            align.alignto(protein, reference, select="backbone")
            results["bb"].append(
                MDAnalysis.analysis.rms.rmsd(
                    reference.select_atoms("backbone").positions,
                    protein.select_atoms("backbone").positions,
                )
            )

        if "composite" in results.keys():
            results["composite"].append(
                results["prmsd"][-1] + (sigma * results["bb"][-1])
            )

    # writing in text files the progress coordinates and the required auxiliary
    # data if needed.
    with open("pcoord.txt", "w") as f:
        for i in range(len(results[settings["pcoord"][0]])):
            line = ""
            for pcoord in settings["pcoord"]:
                line += "{:.4f}    ".format(results[pcoord][i])
            f.write(line + "\n")

    # save fop in file so it can be piped to h5 file
    if "fops" in results.keys():
        if settings["fop_filetype"] == "xyz":
            points_to_xyz_file(
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

    if args.csv is not None:
        import pandas as pd

        results_pd = pd.DataFrame(results)
        results_pd.to_csv(args.csv)
