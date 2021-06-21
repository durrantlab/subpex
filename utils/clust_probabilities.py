"""Summary
This script will read clustering results from cpptraj and obtain final probabilities for each centroid for later use in Boltzmann docking.
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


def check_clustering_matches_h5file(settings, directory):
     # make sure outdir exists if not, creat it
    if not os.path.isdir(args.dir):
        print("Could not load the configuration file with the settings")
        print("make sure the file exists and is correctly formatted")
        raise IOError("Could not load the json file with the settings")
    else:
        pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Script that reads the output of cpptraj clustering and calculates the probabilities of each cluster using the results of the SubPEx run. Outputs a list of probabilities.") 
    parser.add_argument("west", type=str, help="Define the west.h5 file. It is required")
    parser.add_argument("settings", type=str, help="Define the json file with the settings. It is required")
    parser.add_argument("dir", type=str, help="Define the directory were results of docking were stored. It is required")
    args = parser.parse_args()

    # loading west.h5 file
    west = h5py.File(args.west, "r")

    # Loading settings
    settings = check_input_settings_file(args.settings)
    settings = check_clustering_matches_h5file(settings, args.dir)

