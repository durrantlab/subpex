import h5py
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os


def plot_num_walkers(num_walkers, outdir):
    fig = plt.subplots(1, figsize=(10, 7))
    plt.title("Number of Walkers per iteration", fontsize=32)
    plt.plot(range(len(num_walkers)), num_walkers, linewidth=4, color="blue", alpha=0.7)
    plt.tick_params(width=1.0, labelsize=16)
    plt.xlabel("Iteration", fontsize=28)
    plt.ylabel("Number of walkers", fontsize=28)
    plt.savefig(outdir + "/num_walkers.png")


def plot_histogram(results, value, title, directory=None):
    min_1 = None
    max_1 = None
    for i in results[value]:
        if min_1 is None or np.min(i) < min_1:
            min_1 = np.min(i)
        if max_1 is None or np.max(i) > max_1:
            max_1 = np.max(i)

    binsize = list(np.linspace(min_1, max_1, 100))
    fig, axs = plt.subplots(1, 1, figsize=(15, 12))
    axs.set_title("Histogram of {}".format(title), fontsize=26)
    axs.hist(results[value], bins=binsize, histtype="bar", color="blue", alpha=0.7, density=True, stacked=True)
    axs.hist(results[value], bins=binsize, histtype="step", color="black", density=True, stacked=True)
    axs.set_xlabel(title.capitalize(), fontsize=20)
    axs.set_ylabel("PDF", fontsize=20)
    axs.tick_params(width=1.0, labelsize=16)
    if directory is None:
        fig.savefig("{}.png".format(value))
    else:
        fig.savefig(directory + "/{}.png".format(value))


def plot(jd, rmsd, title, outdir, filename):
    fig, axs = plt.subplots(1, 2, figsize=(16, 6))
    fig.suptitle(title, fontsize=32)
    axs[0].set_xlabel("Iteration", fontsize=28)
    axs[1].set_xlabel("Iteration", fontsize=28)
    axs[0].set_ylabel("Jaccard distance", fontsize=28)
    axs[1].set_ylabel("Pocket HA RMSD", fontsize=28)
    axs[0].plot(np.arange(0.5, 0.5 * len(jd) + 0.5, 0.5), jd, label="Jaccard", lw=3.0)
    axs[1].plot(np.arange(0.5, 0.5 * len(jd) + 0.5, 0.5), rmsd, label="RMSD", lw=3.0, color="orange")
    axs[0].tick_params(width=1.0, labelsize=16)
    axs[1].tick_params(width=1.0, labelsize=16)
    plt.savefig(outdir + "/" + filename)


def obtain_paths_higest_rmsd_jd_weight(max_rmsd, max_jd, max_weight, reverse_iterations):
    """

    :param max_rmsd:
    :param max_jd:
    :param max_weight:
    :param reverse_iterations:
    :return:
    """
    rmsd = {}
    rmsd["trajectory"] = [max_rmsd[1]]
    rmsd["jd"] = []
    rmsd["rmsd"] = []
    jd = {}
    jd["trajectory"] = [max_jd[1]]
    jd["jd"] = []
    jd["rmsd"] = []
    weight = {}
    weight["trajectory"] = [max_weight[1]]
    weight["jd"] = []
    weight["rmsd"] = []

    for iteration in reverse_iterations:
        # obtain progress coordinate and parent walker for walker with highest rmsd
        rmsd["jd"].append(west["iterations"][iteration]['pcoord'][rmsd["trajectory"][-1]][-1][0])
        rmsd["rmsd"].append(west["iterations"][iteration]['pcoord'][rmsd["trajectory"][-1]][-1][1])
        rmsd["trajectory"].append(west["iterations"][iteration]['seg_index'][rmsd["trajectory"][-1]][1])
        # obtain progress coordinate and parent walker for walker with highest jaccard distance
        jd["jd"].append(west["iterations"][iteration]['pcoord'][jd["trajectory"][-1]][-1][0])
        jd["rmsd"].append(west["iterations"][iteration]['pcoord'][jd["trajectory"][-1]][-1][1])
        jd["trajectory"].append(west["iterations"][iteration]['seg_index'][jd["trajectory"][-1]][1])
        # obtain progress coordinate and parent walker for walker with highest weight
        weight["jd"].append(west["iterations"][iteration]['pcoord'][weight["trajectory"][-1]][-1][0])
        weight["rmsd"].append(west["iterations"][iteration]['pcoord'][weight["trajectory"][-1]][-1][1])
        weight["trajectory"].append(west["iterations"][iteration]['seg_index'][weight["trajectory"][-1]][1])

    results = {}
    results["rmsd"] = rmsd
    results["jaccard"] = jd
    results["weight"] = weight
    for i in results:
        results[i]["jd"].reverse()
        results[i]["rmsd"].reverse()

    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Script that will read the h5 file and create plots to visualize the progress of the simulation.")
    parser.add_argument("h5file", type=str, help="Define the h5 file file with the settings. It is required")
    parser.add_argument("outdir", type=str, help="Define the output directory to put the files in. It is required")
    parser.add_argument("settings", type=str, help="Define the json file with the settings. It is required")
    args = parser.parse_args()

    try:
        with open(args.settings, "r") as f:
            settings = json.load(f)
    except IOError:
        print("Could not load the json file with the settings")
        print("make sure the file exists and is correctly formatted")
        raise IOError("Could not load the json file with the settings")

    # read h5 file and create outdir
    west = h5py.File(args.h5file, "r")
    os.system("mkdir {}".format(args.outdir))

    # plot num of walkers vs generation
    num_walkers = [len(west["iterations"][k]["pcoord"][()]) for k in west["iterations"].keys()]
    plot_num_walkers(num_walkers, args.outdir)

    reverse_iterations = []
    for i in range(1, len(list(west["iterations"].keys())) + 1):
        reverse_iterations.append(list(west["iterations"].keys())[-i])

    # Obtain walkers with the highest RMSD of pocket atoms and Jaccard distance in the last iteration
    max_rmsd = [0, 0]
    max_jd = [0, 0]
    max_weight = [0, 0]
    for i, value in enumerate(west["iterations"][reverse_iterations[0]]["pcoord"]):
        if value[1][0] > max_jd[0]:
            max_jd = [value[1][0], i]
        if value[1][1] > max_rmsd[0]:
            max_rmsd = [value[1][1], i]
    for i, value in enumerate(west["iterations"][reverse_iterations[0]]['seg_index']):
        if value[0] > max_weight[0]:
            max_weight = [value[0], i]

    results = obtain_paths_higest_rmsd_jd_weight(max_rmsd, max_jd, max_weight, reverse_iterations)
    plot(results["rmsd"]["jd"], results["rmsd"]["rmsd"], "Progress of walker with highest pocket HA RMSD",
         args.outdir, "highest_rmsd.png")
    plot(results["jaccard"]["jd"], results["jaccard"]["rmsd"], "Progress of walker with highest Jaccard distance",
         args.outdir, "highest_jaccard.png")
    plot(results["weight"]["jd"], results["weight"]["rmsd"], "Progress of walker with highest Jaccard distance",
         args.outdir, "highest_weight.png")







