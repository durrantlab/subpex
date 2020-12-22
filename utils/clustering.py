import MDAnalysis as mda
from jdistance import check_input_settings_file
import glob
import argparse
import h5py
import os
import sys
import logging
import numpy as np


def check_clustering_parameters(settings):
    if "clustering" not in settings or type(settings["clustering"]) is not dict:
        logging.critical("There is a problem with the clustering parameters. Please check the settings file.")
        sys.exit("There is a problem with the clustering parameters. Please check the settings file.")
    else:
        pass
    # todo finish the function


def get_bins_dictionary(west, settings, max_clusters_per_bin, final_clusters):
    """
    get_bins_dictionary is a function that created a dictionary with unique keys for each occupied bin. And the elements
    are the walker generation/iteration and the walker number.
    :param west: The h5 file containing the results of the SubPEx (weighted ensemble) simulation.
    :param settings: the settings dictionary.
    :param max_clusters_per_bin: maximum number of clusters to calculate per bin.
    :param final_clusters: number of final clusters.
    :return bins: A dictionary containing walker identifiers per bin.
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


def concatenate_trajectories(bin_dictionary, directory, num_clust):
    pass


def get_clustering_command_cpptraj(method, num_clusters, selection_string):
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


def cluster_cpptaj(bins, settings, dir, residues_pocket):
    """

    :param bins:
    :param settings:
    :param dir:
    :param residues_pocket:
    :return:
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

def make_cpptraj_pdb_extractor(walkers):
    """
    This is a function that ...
    :param walkers:
    :return:
    """
    text = "blah"
    return text


def get_selection_cpptraj(residues_pocket):
    """
    get_selection_cpptraj is a function that takes a list of residue numbers and returns a selection string formated for
    cpptraj. The selection will take all heavy atoms of each of the residues.
    :param residues_pocket: list of residue numbers
    :return: selection string formatted for cpptraj
    """
    selection_string = "'("
    for resid in residues_pocket[:-1]:
        selection_string += ":{}|".format(resid)
    selection_string += ":{})&!@H'".format(residues_pocket[-1])
    return selection_string


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Script that reads the concateneted dcd files and performs clustering.")
    parser.add_argument("west", type=str, help="Define the west.h5 file. It is required")
    parser.add_argument("settings", type=str, help="Define the json file with the settings. It is required")
    parser.add_argument("dir", type=str, help="Define the directory were results will be stored. It is required")
    parser.add_argument("--num_clust", type=int, default=25, help="Define the final number of clusters. Default is 25")
    parser.add_argument("--max_clust_num", type=int, default=25, help="Define the maximum number of clusters to obtain per bin. Default is 25.")
    parser.add_argument("--num_proc", type=int, default=-1, help="Define the number of processors to use. Default max amount - 2.")
    parser.add_argument("--cpptraj", action="store_true", help="If you want the input files created for CPPTRAJ.")
    args = parser.parse_args()

    # loading west.h5 file
    west = h5py.File(args.west, "r")

    # Loading settings
    settings = check_input_settings_file(args.settings)

    # Make sure the directory we will store files exists
    directory = glob.glob(args.dir)
    if len(directory) == 0:
        os.system("mkdir {}".format(args.dir))
    else:
        pass

    # Get residues that form the pocket
    with open(settings["selection_file"], "r") as f:
        selection_pocket = f.readlines()[0]
    text_split = selection_pocket.split("or")
    text_split[len(text_split) - 1] = text_split[-1].strip("and not name H") # todo change this part
    residues_pocket = [x.strip(" ") for x in text_split]
    residues_pocket = [int(x.strip("resid")) for x in residues_pocket]

    # sdjjfaskdjf
    bins = get_bins_dictionary(west, settings, args.max_clust_num, args.num_clust)

    if args.cpptraj:
        cluster_cpptaj(bins, settings, args.dir, residues_pocket)
    else:
        cluster_per_bin(bins, dir, settings["west_home"])

    #concatenate_trajectories(bins, args.dir,)
