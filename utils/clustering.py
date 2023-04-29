"""Summary
This script will generate files to cluster a SubPEx run
"""

from typing import List
import MDAnalysis as mda
import glob
import argparse
import contextlib
import h5py
import os
import sys
import logging
import numpy as np

currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from westpa_scripts.pcoord import check_input_settings_file


def check_clustering_parameters(settings: dict) -> dict:
    """Checks the clustering parameters in the settings.

    Args:
        settings (dict): The settings dictionary.

    Returns:
        dict: The settings dictionary.
    """

    if "clustering" not in settings or type(settings["clustering"]) is not dict:
        logging.critical(
            "There is a problem with the clustering parameters. Please check the settings file."
        )
        sys.exit(
            "There is a problem with the clustering parameters. Please check the settings file."
        )

    # TODO finish check_clustering_parameters
    return settings


def get_bins_dictionary(west: dict, settings: dict) -> dict:
    """Creates a dictionary with the bins and the walkers in each bin.
    
    Args:
        west (dict): The west.h5 file as a dictionary.
        settings (dict): The settings dictionary.
        
    Returns:
        dict: A dictionary with the bins and the walkers in each bin.
    """

    # create an empty dictionary to store the bins info
    bins = {}
    # now place each walker into a bin
    for iteration in west["iterations"]:
        for walker in range(len(west["iterations"][iteration]["pcoord"])):
            name = ""
            for key in settings["clustering"]["bins"]:
                if key in settings["pcoord"]:
                    for i in range(len(settings["pcoord"])):
                        if key == settings["pcoord"][i]:
                            value = west["iterations"][iteration]["pcoord"][-1][-1][i]
                elif key in settings["auxdata"]:
                    value = west["iterations"][iteration]["auxdata"][key][walker]
                for a, bin_val in enumerate(settings["clustering"]["bins"][key]):
                    if value <= bin_val:
                        name += f"{a}_"
                        break
            name = name[:-1]  # I only need to take the last underscore
            if name in bins:
                bins[name]["walkers"].append((iteration, walker))
            else:
                bins[name] = {"walkers": [(iteration, walker)]}
    return bins


def calculate_clusters_per_bin(bins: dict, settings: dict):
    """Calculates the number of clusters per bin. Updates bin in place.

    Args:
        bins (dict): A dictionary with the bins and the walkers in each bin.
        settings (dict): The settings dictionary.
    """

    # determine the maximum number of walkers in a bin for normalization of clusters per bin
    max_num_walkers = 0
    for key in bins:
        if len(bins[key]["walkers"]) > max_num_walkers:
            max_num_walkers = len(bins[key])
    # append as the last element of the list the number of clusters we will obtain in said bin
    for key in bins:
        if (
            len(bins[key]["walkers"])
            <= settings["clustering"]["min_number_clusters_generation_bin"]
        ):
            bins[key]["clusters"] = len(bins[key]["walkers"])
        elif (len(bins[key]["walkers"]) / max_num_walkers) * settings["clustering"][
            "max_number_clusters_generation_bin"
        ] <= 3:
            bins[key]["clusters"] = 3
        else:
            value = int(
                round(
                    (len(bins[key]["walkers"]) / max_num_walkers)
                    * settings["clustering"]["max_number_clusters_generation_bin"],
                    0,
                )
            )
            bins[key]["clusters"] = value


def create_bin_cpptraj_files(bins: dict, settings: dict, directory: str):
    """Creates the cpptraj files for each bin.

    Args:
        bins (dict): A dictionary with the bins and the walkers in each bin.
        settings (dict): The settings dictionary.
        directory (str): The directory where the cpptraj files will be created.
    """

    calculate_clusters_per_bin(bins, settings)
    # will generate the selection string that will be used
    if settings["clustering"]["clustering_region"] == "pocket":
        selection_string = get_selection_cpptraj(settings["selection_file"])
    elif settings["clustering"]["clustering_region"] == "backbone":
        selection_string = "'@CA,C,O,N'"

    for key in bins:
        os.system(f"mkdir {directory}/{key}")
        create_cpptraj_file_bins(bins, settings, directory, key, selection_string)


def create_cpptraj_file_bins(
    bins: dict, settings: dict, directory: str, key: str, selection_string: str
):
    """Creates the cpptraj file for a bin.

    Args:
        bins (dict): A dictionary with the bins and the walkers in each bin.
        settings (dict): The settings dictionary.
        directory (str): The directory where the cpptraj files will be created.
        key (str): The name of the bin.
        selection_string (str): The selection string that will be used.
    """

    with open(f"{directory}/{key}/clustering.in", "w") as f:
        f.write(f'parm {settings["topology"]} \n')
        for walker in bins[key]["walkers"]:
            filename = glob.glob(settings["west_home"])[0]
            iteration = int(walker[0].split("_")[1])
            walker_num = walker[1]
            filename += "traj_segs/{:06}/{:06}/seg.dcd".format(iteration, walker_num)
            if os.path.exists(filename):
                f.write(f"trajin {filename} \n")
        clustering_command = get_clustering_command_cpptraj(
            settings["clustering"], bins[key]["walkers"]["clusters"], selection_string
        )
        f.write(clustering_command)
        f.write("\n")


def get_selection_cpptraj(filename: str) -> str:
    """Takes the filename of the selection string used in SubPEx and returns a
    selection string formated for cpptraj. The selection will take all heavy
    atoms of each of the residues.

    Args:
        filename (str): The filename of the selection string used in SubPEx.

    Returns:
        str: The selection string formated for cpptraj.
    """

    # open selection files used in subpex
    with open(filename, "r") as f:
        selection_mdanalysis = f.readlines()[0]

    residues_pocket = []

    for i in selection_mdanalysis.split():
        with contextlib.suppress(Exception):
            residues_pocket.append(int(i))

    # Sorting it out just because I wanted
    residues_pocket.sort()

    # converting the residues list to selection screen used in cpptraj
    selection_string = "'("
    for resid in residues_pocket[:-1]:
        selection_string += f":{resid}|"
    selection_string += f":{residues_pocket[-1]})&!@H'"
    return selection_string


def get_clustering_generation_cpptraj(west_file: dict, settings: dict, directory: str):
    """Creates the cpptraj files for each bin.

    Args:
        west_file (dict): The west.h5 file contents.
        settings (dict): The settings dictionary.
        directory (str): The directory where the cpptraj files will be created.
    """

    # will generate the selection string that will be used
    if settings["clustering"]["clustering_region"] == "pocket":
        selection_string = get_selection_cpptraj(settings["selection_file"])
    elif settings["clustering"]["clustering_region"] == "backbone":
        selection_string = "'@CA,C,O,N'"

    west_home = glob.glob(settings["west_home"])[0]
    all_in_files = []
    number_clusters = get_number_clusters_generation(
        west_file, settings["clustering"]["max_number_clusters_generation_bin"]
    )
    if not os.path.isdir(f"{directory}/last_clustering"):
        os.system(f"mkdir {directory}/last_clustering")
    for i, iteration in enumerate(west_file["iterations"].keys()):
        directory_gen = f'{west_home}traj_segs/{iteration.split("_")[1][-6:]}'
        generation_number = directory_gen.split("/")[-1]
        # make sure the directory exists (this is because the last iter that appears on the
        # h5 file has not yet started)
        if os.path.isdir(directory_gen):
            walkers = get_generation_walkers_list(directory_gen)
        else:
            continue
        # now create a directory where the input file for cpptraj will be created.
        if not os.path.isdir(directory + generation_number):
            os.system(f"mkdir {directory}/{generation_number}")
        all_in_files.append(generation_number)
        make_input_file_cpptraj(
            f"{directory}/{generation_number}",
            walkers,
            settings,
            number_clusters[i],
            selection_string,
        )
    with open(f"{directory}/run_cpptraj_per_generation.sh", "w") as f:
        for i in all_in_files:
            f.write(f"cd {i} \n")
            f.write("cpptraj.OMP -i clustering.in > clustering.log \n")
            f.write("cd ../ \n")
        f.write("cd last_clustering \n")
        f.write("cpptraj.OMP -i clustering.in > clustering.log \n")
    with open(f"{directory}/last_clustering/clustering.in", "w") as f:
        f.write(f'parm {settings["topology"]} \n')
        for i, numgen in enumerate(all_in_files):
            for j in range(number_clusters[i]):
                f.write(f"trajin ../{numgen}/rep.c{j}.pdb \n")
        text = get_clustering_command_cpptraj(
            settings["clustering"]["method"],
            settings["clustering"]["number_clusters"],
            selection_string,
        )
        for i in text:
            f.write(i)


def get_generation_walkers_list(directory_gen: str) -> List[str]:
    """Gets the list of walkers for a given generation.

    Args:
        directory_gen (str): The directory where the generation is.

    Returns:
        List[str]: The list of walkers.
    """

    # will check later what the extension for coordinate files is. will work for dcd and nc.
    file_extension = None
    walkers_list = []
    # now search for all the seg.dcd or seg.nc files and append them to walker_lsit
    walkers = glob.glob(f"{directory_gen}/*/")
    for walker in walkers:
        if file_extension is None:
            if os.path.isfile(f"{walker}/seg.dcd"):
                file_extension = "dcd"
            elif os.path.isfile(f"{walker}/seg.nc"):
                file_extension = "nc"
        if os.path.isfile(f"{walker}seg.{file_extension}"):
            walkers_list.append(f"{walker}seg.{file_extension}")
    walkers_list.sort()
    return walkers_list


def make_input_file_cpptraj(
    directory: str,
    walkers_list: List[str],
    settings: dict,
    num_clusters: int,
    selection_string: str,
):
    """Creates the input file for cpptraj.

    Args:
        directory (str): The directory where the input file will be created.
        walkers_list (List[str]): The list of walkers for the generation.
        settings (dict): The settings dictionary.
        selection_string (str): The selection string for cpptraj.
    """

    text = f'parm {settings["topology"]} \n'
    for walker in walkers_list:
        text += f"trajin {walker} \n"
    text += get_clustering_command_cpptraj(
        settings["clustering"]["method"], num_clusters, selection_string
    )
    with open(f"{directory}/clustering.in", "w") as f:
        f.write(text)


def get_clustering_command_cpptraj(
    method: str, num_clusters: int, selection_string: str
) -> str:
    """Creates the clustering command for cpptraj.
    
    Args:
        method (str): the clustering method to use. Currently, only hierarchical is implemented.
        num_clusters (int): The number of clusters to use.
        selection_string (str): The selection string for cpptraj.
    
    Returns:
        str: The clustering command for cpptraj.
    """

    if method == "hierarchical":
        text = """cluster C0 \ 
    hieragglo clusters {num_clust} averagelinkage epsilonplot epsilonplot.dat \ 
    dme {selection} \ 
    out cpptraj_cluster.dat \ 
    savepairdist pairdist pairdist.cpp \ 
    sil Sil \ 
    summary summary.dat \ 
    info info.dat \ 
    cpopvtime cpopvtime.agr normframe \ 
    repout rep repfmt pdb singlerepout singlerep.nc singlerepfmt netcdf \ 
    avgout Avg avgfmt restart""".format(
            num_clust=num_clusters, selection=selection_string
        )
    else:
        print("Have not implemented this method yet")  # todo add other methods
        sys.exit(
            "Please change clustering algorithm or adapt script to accommodate the clustering method"
        )
    return text


def get_number_clusters_generation(
    west_file: dict, max_clusters: int, min_clusters: int = 3
) -> List[int]:
    """Gets the number of clusters for each generation.

    Args:
        west_file (dict): The west.h5 file.
        max_clusters (int): The maximum number of clusters to use.
        min_clusters (int, optional): The minimum number of clusters to use.
            Defaults to 3.

    Returns:
        List[int]: The number of clusters for each generation.
    """

    number_walkers = [
        len(west["iterations"][iteration]["seg_index"])
        for iteration in west_file["iterations"].keys()
    ]
    if min_clusters > np.min(number_walkers):
        min_clusters = np.min(number_walkers)
    max_walkers = np.max(number_walkers)
    number_clusters = []
    for i in number_walkers:
        val = int(np.round((i / max_walkers) * 25))
        val = max(val, min_clusters)
        number_clusters.append(val)
    return number_clusters


def get_clustering_bins_cpptraj(west: dict, settings: dict, directory: str):
    """Creates the input files (sh) for cpptraj to cluster the bins.

    Args:
        west (dict): The west.h5 file.
        settings (dict): The settings dictionary.
        directory (str): The directory where the input files will be created.
    """

    bins = get_bins_dictionary(west, settings)
    create_bin_cpptraj_files(bins, settings, directory)
    with open(f"{directory}/run_clustering_bins.sh", "w") as f:
        for key in bins.keys:
            f.write(f"cd {directory}/{key} \n")
            f.write("cpptraj.OMP -i clustering.in > clustering.log \n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Prepares files to run clustering on a SubPEx run via cpptraj. See 'clustering' section of west.cfg to configure further (e.g., cluster per generation or dividing the space in bins)."
    )
    parser.add_argument(
        "west", type=str, help="Path to the west.h5 file. It is required"
    )
    parser.add_argument(
        "settings",
        type=str,
        help="Path to the yaml file with the settings (e.g., west.cfg). It is required",
    )
    parser.add_argument(
        "dir",
        type=str,
        help="Path to the directory where files will be saved (e.g., ./cluster/). It is required",
    )
    args = parser.parse_args()

    # loading west.h5 file
    west = h5py.File(args.west, "r")

    # Loading settings
    settings = check_input_settings_file(args.settings)
    settings = check_clustering_parameters(settings)

    # make sure outdir exists if not, creat it
    if not os.path.isdir(args.dir):
        os.system(f"mkdir {args.dir}")

    # Obtaining full path to directory
    directory = glob.glob(args.dir)[0]

    if settings["clustering"]["clustering_engine"] == "cpptraj":
        if settings["clustering"]["cluster_generation"]:
            get_clustering_generation_cpptraj(west, settings, directory)
        elif settings["clustering"]["cluster_bin"]:
            get_clustering_bins_cpptraj(west, settings, directory)

    elif settings["clustering"]["clustering_engine"] == "MDAnalysis":
        print("Not yet implemented")

    print("\nDone! See files in " + directory)
