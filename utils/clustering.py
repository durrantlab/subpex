"""Summary
This script will generate files for clustering or cluster itself the output of a SubPEx run
"""
import MDAnalysis as mda
#from .jdistance import check_input_settings_file
import glob
import argparse
import h5py
import os
import sys
import logging
import numpy as np
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from westpa_scripts.jdistance import check_input_settings_file

def check_clustering_parameters(settings):
    """[summary]

    Args:
        settings ([type]): [description]
    """
    if "clustering" not in settings or type(settings["clustering"]) is not dict:
        logging.critical("There is a problem with the clustering parameters. Please check the settings file.")
        sys.exit("There is a problem with the clustering parameters. Please check the settings file.")
    else:
        pass
    #TODO finish check_clustering_parameters
    return settings


def get_bins_dictionary(west, settings, max_clusters_per_bin, final_clusters):
    """
    get_bins_dictionary is a function that created a dictionary with unique keys for each occupied bin. And the elements
    are the walker generation/iteration and the walker number.
    :param west: The h5 file containing the results of the SubPEx (weighted ensemble) simulation.
    :param settings: the settings dictionary.
    :param max_clusters_per_bin: maximum number of clusters to calculate per bin.
    :param final_clusters: number of final clusters.
    :return bins: A dictionary containing walker identifiers per bin.
    
    Args:
        west (TYPE): Description
        settings (TYPE): Description
        max_clusters_per_bin (TYPE): Description
        final_clusters (TYPE): Description
    
    Returns:
        TYPE: Description
    """
    # initialize an empty directory which will contain all the occupied bins and walker per bin.
    bins = {}
    # iterate over all generations and walkers to save them into the bin dictionary
    for iteration in west["iterations"]:
        for walker, pcoord in enumerate(west["iterations"][iteration]["pcoord"]):
            # checking Jaccard distance
            for i, value in enumerate(settings["clustering"]["bins"]["jd"]):
                if pcoord[-1][0] < value:
                    break
            # checking pocket HA RMSD
            for n, value in enumerate(settings["clustering"]["bins"]["rmsd"]):
                if pcoord[-1][1] < value:
                    break
            # check if it is the first element in the bin
            if "bin_{}_{}".format(i, n) in bins:
                bins["bin_{}_{}".format(i, n)].append([iteration, walker])
            else:
                bins["bin_{}_{}".format(i, n)] = []
                bins["bin_{}_{}".format(i, n)].append([iteration, walker])

    # calculate how many clusters we will get for each bin (it will be proportional to the number of walkers)
    calculate_clusters_per_bin(bins, max_clusters_per_bin)
    bins["final_clusters"] = final_clusters

    return bins


def calculate_clusters_per_bin(bins, max_clusters_per_bin, minimum_number_clusters=3):
    """
    calculate_clusters_per_bin is a function that takes teh bins dictionary, the maximum number of clusters to obtain
    per bin and can take the minimum number of clusters per bin and appends to each bin the number of clusters to
    calculate for that specific bin.
    :param bins: dictionary with bins list containing walker generation and walker number in said bin.
    :param max_clusters_per_bin: int maximum number of clusters to obtain per bin,
    :param minimum_number_clusters: int the minimum number of clusters to calculate per bin.
    :return:
    
    Args:
        bins (TYPE): Description
        max_clusters_per_bin (TYPE): Description
        minimum_number_clusters (int, optional): Description
    """
    # determine the maximum number of walkers in a bin for normalization of clusters per bin
    max_num_walkers = 0
    for key in bins.keys():
        if len(bins[key]) > max_num_walkers:
            max_num_walkers = len(bins[key])
        else:
            pass

    # append as the last element of the list the number of clusters we will obtain in said bin
    for key in bins.keys():
        if len(bins[key]) <= minimum_number_clusters:
            bins[key].append(len(bins[key]))
        elif (len(bins[key]) / max_num_walkers) * max_clusters_per_bin <= 3:
            bins[key].append(3)
        else:
            value = int(round((len(bins[key]) / max_num_walkers) * max_clusters_per_bin, 0))
            bins[key].append(value)


def cluster_cpptaj(bins, settings, dir, residues_pocket):
    """
    
    :param bins:
    :param settings:
    :param dir:
    :param residues_pocket:
    :return:
    
    Args:
        bins (TYPE): Description
        settings (TYPE): Description
        dir (TYPE): Description
        residues_pocket (TYPE): Description
    """
    # Create text for final clustering script, and list containing all the files so it can be written later.
    final_clustering = "parm {} \n".format(settings["topology"])
    all_in_files = []
    # obtain selection string for cpptraj from residue numbers
    selection_string = get_selection_cpptraj(residues_pocket)
    # loop over bins to create a the directory and infile for clustering
    for bin in bins.keys():
        # make directory for the bin, results of clustering will be saved here.
        os.system("mkdir {}/{}".format(dir, bin))
        # avoid final cluster since it is only a holder for number of clusters
        if bin != "final_clusters":
            # add file to all_in_files for bash script
            all_in_files.append("{}/{}/cluster_{}.in".format(dir, bin, bin))
            # crate input file for clustering in a specific bin
            with open("{}/{}/cluster_{}.in".format(dir, bin, bin), "w") as f:
                # load parameter file
                f.write("parm {} \n".format(settings["topology"]))
                # the last element of the list is the number of clusters to calculate
                num_clusters = bins[bin].pop()
                # make sure we have more walkers than the necessary number of clusters
                if num_clusters >= len(bins[bin]):
                    # if not, just extract the last frame from each trajectory
                    text = make_cpptraj_pdb_extractor(bins[bin])
                else:
                    text = ""
                    # obtain the text for the clustering commands for the cpptraj infile.
                    for walker in bins[bin]:
                        iteration = walker[0].split("_")[1]
                        iteration = iteration[np.absolute(len(iteration) - 6):]
                        walker_number = (np.absolute(6 - len(str(walker[1])))) * "0" + str(walker[1])
                        filename = "{}/{}/{}/{}/seg.dcd".format(settings["west_home"], "traj_segs", iteration,
                                                                       walker_number)
                        # make sure the file exists (may be redundant)
                        if len(glob.glob(filename)) == 1:
                            text += ("trajin {} \n".format(filename))
                    text += get_clustering_command_cpptraj(settings["clustering"]["method"], num_clusters,
                                                          selection_string)
                f.write(text)
                bins[bin].append(num_clusters)
            for i in range(num_clusters):
                final_clustering += "trajin {}/{}/rep.c{}.pdb \n".format(dir, bin, i)
    text = get_clustering_command_cpptraj(settings["clustering"]["method"], bins["final_clusters"], selection_string)
    final_clustering += text
    with open("run_cpptraj_per_bin.sh", "w") as f:
        for i in all_in_files:
            f.write("cpptraj.OMP -i {} > {}.log \n".format(i, i.split(".")[0]))
    os.system("mkdir {}/final_clustering".format(dir))
    with open("{}/final_clustering/clustering.in".format(dir), "w") as f:
        f.write(final_clustering)


def get_selection_cpptraj(filename):
    """
    get_selection_cpptraj is a function that takes the filename of the selection string used in 
    SubPEx and returns a selection string formated for cpptraj. The selection will take all 
    heavy atoms of each of the residues.

    Args:
        filename ([type]): [description]

    Returns:
        [type]: [description]
    """
    # open selection files used in subpex
    with open(filename, "r") as f:
        selection_mdanalysis = f.readlines()[0]

    residues_pocket = []

    for i in selection_mdanalysis.split():
        try:
            residues_pocket.append(int(i))
        except:
            pass
    
    # Sorting it out just because I wanted
    residues_pocket.sort()

    # converting the residues list to selection screen used in cpptraj
    selection_string = "'("
    for resid in residues_pocket[:-1]:
        selection_string += ":{}|".format(resid)
    selection_string += ":{})&!@H'".format(residues_pocket[-1])
    return selection_string


def get_clustering_generation_cpptraj(west_file, settings, directory):
    # will generate the selection string that will be used
    if settings["clustering"]["clustering_region"] == "pocket":
        selection_string = get_selection_cpptraj(settings["selection_file"])
    elif settings["clustering"]["clustering_region"] == "backbone":
        selection_string = "'@CA,C,O,N'"
    else:
        pass
    west_home = glob.glob(settings["west_home"])[0]
    all_in_files = []
    number_clusters = get_number_clusters(west_file, settings["clustering"]["max_number_clusters_generation_bin"])
    if not os.path.isdir(directory + "/last_clustering"):
        os.system("mkdir {}".format(directory + "/last_clustering" ))
    for i, iteration in enumerate(west_file["iterations"].keys()):
        directory_gen = west_home + "traj_segs/{}".format(iteration.split("_")[1][-6:])
        generation_number = directory_gen.split("/")[-1]
        # make sure the directory exists (this is because the last iter that appears on the 
        # h5 file has not yet started)
        if os.path.isdir(directory_gen):
            walkers = get_generation_walkers_list(directory_gen)
        else:
            continue
        # now create a directory where the input file for cpptraj will be created.
        if not os.path.isdir(directory + generation_number):
            os.system("mkdir {}".format(directory + "/" + generation_number))
        all_in_files.append(generation_number)
        make_input_file_cpptraj(directory + "/" + generation_number, walkers, settings, number_clusters[i], selection_string)
    with open(directory + "/run_cpptraj_per_generation.sh", "w") as f:
        for i in all_in_files:
            f.write("cd {} \n".format(i)) 
            f.write("cpptraj.OMP -i clustering.in > clustering.log \n")
            f.write("cd ../ \n")
        f.write("cd last_clustering \n")
        f.write("cpptraj.OMP -i last_clustering > last_clustering.log \n")
    with open(directory + "/last_clustering/clustering.in", "w") as f:
        f.write("parm {} \n".format(settings["topology"]))
        for i, numgen in enumerate(all_in_files):
            for j in range(number_clusters[i]):
                f.write("trajin ../{}/rep.c{}.pdb \n".format(numgen, j))
        text = get_clustering_command_cpptraj(settings["clustering"]["method"], 
                                   settings["clustering"]["number_clusters"], selection_string)
        for i in text:
            f.write(i)
            

def get_generation_walkers_list(directory_gen):
    # will check later what the extension for coordinate files is. will work for dcd and nc.
    file_extension = None
    walkers_list = []
    # now search for all the seg.dcd or seg.nc files and append them to walker_lsit
    walkers = glob.glob(directory_gen + "/*/")
    for walker in walkers:
        if file_extension is None:
            if os.path.isfile(walker + "/seg.dcd"):
                file_extension = "dcd"
            elif os.path.isfile(walker + "/seg.nc"):
                file_extension = "nc"
        if os.path.isfile(walker + "seg.{}".format(file_extension)):
            walkers_list.append(walker + "seg.{}".format(file_extension))
    walkers_list.sort()
    return walkers_list


def make_input_file_cpptraj(directory, walkers_list, settings, num_clusters, selection_string):
    """[summary]

    Args:
        directory ([type]): [description]
        walkers_list ([type]): [description]
        settings ([type]): [description]
        selection_string ([type]): [description]
    """
    text = "parm {} \n".format(settings["topology"])
    for walker in walkers_list:
        text += "trajin {} \n".format(walker)
    text += get_clustering_command_cpptraj(settings["clustering"]["method"], num_clusters, selection_string)
    with open(directory + "/clustering.in", "w") as f:
        f.write(text)


def get_clustering_command_cpptraj(method, num_clusters, selection_string):
    """Summary
    
    Args:
        method (TYPE): Description
        num_clusters (TYPE): Description
        selection_string (TYPE): Description
    
    Returns:
        TYPE: Description
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
    avgout Avg avgfmt restart""".format(num_clust=num_clusters, selection=selection_string)
    else:
        print("Have not implemented this method yet")  # todo add other methods
        sys.exit("Please change clustering algorithm or adapt script to accommodate the clustering method")
    return text


def get_number_clusters(west_file, max_clusters, min_clusters=3):
    number_walkers = []
    for iteration in west_file["iterations"].keys():
        number_walkers.append(len(west["iterations"][iteration]['seg_index']))
    if min_clusters > np.min(number_walkers):
        min_clusters = np.min(number_walkers)    
    max_walkers = np.max(number_walkers)
    number_clusters = []
    for i in number_walkers:
        val = int(np.round((i / max_walkers) * 25))
        if val < min_clusters:
            val = min_clusters
        number_clusters.append(val)
    return number_clusters


def get_clustering_bins_cpptraj(west, settings):
    print("hiiiii")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Script that reads the concateneted dcd files and performs clustering.")
    parser.add_argument("west", type=str, help="Define the west.h5 file. It is required")
    parser.add_argument("settings", type=str, help="Define the json file with the settings. It is required")
    parser.add_argument("dir", type=str, help="Define the directory were results will be stored. It is required")
    args = parser.parse_args()

    # loading west.h5 file
    west = h5py.File(args.west, "r")

    # Loading settings
    settings = check_input_settings_file(args.settings)
    settings = check_clustering_parameters(settings)

    # make sure outdir exists if not, creat it
    if not os.path.isdir(args.dir):
        os.system("mkdir {}".format(args.dir))
    
    # Obtaining full path to directory
    directory = glob.glob(args.dir)[0] 

    if settings["clustering"]["clustering_engine"] == "cpptraj":
        if settings["clustering"]["cluster_generation"]:
            get_clustering_generation_cpptraj(west, settings, directory)
        elif settings["clustering"]["cluster_bin"]:
            get_clustering_bins_cpptraj(west, settings, directory)
        else:
            pass
    elif settings["clustering"]["clustering_engine"] == "MDAnalysis":
        print("Not yet implemented")
