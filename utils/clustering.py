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
    if "clustering" not in settings["subpex"] or type(settings["subpex"]["clustering"]) is not dict:
        logging.critical("There is a problem with the clustering parameters. Please check the settings file.")
        sys.exit("There is a problem with the clustering parameters. Please check the settings file.")
    else:
        pass
    # TODO finish check_clustering_parameters
    return settings


def get_bins_dictionary(west, settings):
    # create an empty dictionary to store the bins info
    bins = {}
    # now place each walker into a bin
    for iteration in west["iterations"]:
        for walker in range(len(west["iterations"][iteration]["pcoord"])):
            name = ""
            for key in settings["subpex"]["clustering"]["bins"]:
                if key in settings["subpex"]["pcoord"]:
                    for i in range(len(settings["subpex"]["pcoord"])):
                        if key == settings["subpex"]["pcoord"][i]:
                            value = west["iterations"][iteration]["pcoord"][-1][-1][i]
                elif key in settings["subpex"]["auxdata"]:
                    value = west["iterations"][iteration]["auxdata"][key][walker]
                for a, bin_val in enumerate(settings["subpex"]["clustering"]["bins"][key]):
                    if value > bin_val:
                        pass
                    else:
                        name += "{}_".format(a)
                        break
            name  = name[:-1] # I only need to take the last underscore
            if name in bins.keys():
                bins[name]["walkers"].append((iteration, walker))
            else:
                bins[name] = {"walkers": [(iteration, walker)]}
    return bins
    

def calculate_clusters_per_bin(bins, settings):
    # determine the maximum number of walkers in a bin for normalization of clusters per bin
    max_num_walkers = 0
    for key in bins.keys():
        if len(bins[key]["walkers"]) > max_num_walkers:
            max_num_walkers = len(bins[key])
        else:
            pass
    # append as the last element of the list the number of clusters we will obtain in said bin
    for key in bins.keys():
        if len(bins[key]["walkers"]) <= settings["subpex"]["clustering"]["min_number_clusters_generation_bin"]:
            bins[key]["clusters"] = len(bins[key]["walkers"])
        elif (len(bins[key]["walkers"]) / max_num_walkers) * settings["subpex"]["clustering"]["max_number_clusters_generation_bin"] <= 3:
            bins[key]["clusters"] = 3
        else:
            value = int(round((len(bins[key]["walkers"]) / max_num_walkers) * settings["subpex"]["clustering"]["max_number_clusters_generation_bin"], 0))
            bins[key]["clusters"] = value


def create_bin_cpptraj_files(bins, settings, directory):
    calculate_clusters_per_bin(bins, settings)
    # will generate the selection string that will be used
    if settings["clustering"]["clustering_region"] == "pocket":
        selection_string = get_selection_cpptraj(settings["selection_file"])
    elif settings["clustering"]["clustering_region"] == "backbone":
        selection_string = "'@CA,C,O,N'"
    else:
        pass
    for key in bins.keys():
        os.system("mkdir {}".format(directory + "/{}".format(key)))
        create_cpptraj_file_bins(bins, settings, directory, key, selection_string)


def create_cpptraj_file_bins(bins, settings, directory, key, selection_string):
    with open(directory + "/{}/clustering.in".format(key), "w") as f:
        f.write("parm {} \n".format(settings["subpex"]["topology"]))
        for walker in bins[key]["walkers"]:
            filename = glob.glob(settings["west_home"])[0]
            iteration = int(walker[0].split("_")[1])
            walker_num = walker[1]
            filename += "traj_segs/{:06}/{:06}/seg.dcd".format(iteration, walker_num)
            if os.path.exists(filename):
                f.write("trajin {} \n".format(filename))
        clustering_command = get_clustering_command_cpptraj(settings["clustering"], bins[key]["walkers"]["clusters"], selection_string)
        f.write(clustering_command)
        f.write("\n")
        

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
    number_clusters = get_number_clusters_generation(west_file, settings["clustering"]["max_number_clusters_generation_bin"])
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
        f.write("cpptraj.OMP -i clustering.in > clustering.log \n")
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


def get_number_clusters_generation(west_file, max_clusters, min_clusters=3):
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


def get_clustering_bins_cpptraj(west, settings, directory):
    bins = get_bins_dictionary(west, settings)
    create_bin_cpptraj_files(bins, settings, directory)
    for 

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Script that performs clustering or prepares files to run clustering using cpptraj. It can cluster per generation or divided the space in bins according to what the user wants.")
    parser.add_argument("west", type=str, help="Define the west.h5 file. It is required")
    parser.add_argument("settings", type=str, help="Define the yaml file with the settings. It is required")
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
