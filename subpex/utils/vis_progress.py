from typing import List

import argparse
import os

import h5py
import matplotlib.pyplot as plt
import numpy as np

from .westpa_scripts.jdistance import check_file_exists, check_input_settings_file


def rolling_average(data: List[float], window: int) -> List[float]:
    """Gets the rolling average. The end of the list will only make the window
    smaller.

    Args:
        data (list): list of the data to obtain the rolling average.
        window (int): size of the window to use.

    Returns:
        average (list): list containing the rolling average.
    """

    average = []
    for i in range(len(data)):
        if i + window < len(data):
            average.append(np.mean(data[i : i + window]))
        else:
            average.append(np.mean(data[i:]))

    return average


def plot_num_walkers(west: h5py.File, outdir: str) -> None:
    """Plots the number of walkers per iteration. The plot will be saved to the
    specified directory.

    Args:
        west (h5file): h5 file containing all the info for the simulation.
        outdir (str): path to outdir.
    """

    # get number of walkers on each iteration
    num_walkers = [
        len(west["iterations"][k]["pcoord"][()]) for k in west["iterations"].keys()
    ]
    # Plotting the num of walkers per iteration.
    fig = plt.subplots(1, figsize=(10, 7))
    plt.title("Number of Walkers per iteration", fontsize=32)
    plt.plot(
        range(1, len(num_walkers) + 1),
        num_walkers,
        linewidth=3,
        color="blue",
        alpha=0.7,
    )
    plt.tick_params(width=1.0, labelsize=16)
    plt.xlabel("Iteration", fontsize=28)
    plt.ylabel("Number of walkers", fontsize=28)
    plt.savefig(f"{outdir}/num_walkers.png")


def plot_histogram(results: dict, key: str, title: str, directory: str = None) -> None:
    """Plots the histogram of the specified key. The plot will be saved to the
    specified directory.

    Args:
        results (dict): dictionary containing the results of the analysis.
        key (str): key to plot.
        title (str): title of the plot.
        directory (str): path to outdir.
    """

    # TODO: Erich - refactor plot_histogram function
    min_1 = None
    max_1 = None
    for i in results[key]:
        if min_1 is None or np.min(i) < min_1:
            min_1 = np.min(i)
        if max_1 is None or np.max(i) > max_1:
            max_1 = np.max(i)

    binsize = list(np.linspace(min_1, max_1, 100))
    fig, axs = plt.subplots(1, 1, figsize=(15, 12))
    axs.set_title(f"Histogram of {title}", fontsize=26)
    axs.hist(
        results[key],
        bins=binsize,
        histtype="bar",
        color="blue",
        alpha=0.7,
        density=True,
        stacked=True,
    )
    axs.hist(
        results[key],
        bins=binsize,
        histtype="step",
        color="black",
        density=True,
        stacked=True,
    )
    axs.set_xlabel(title.capitalize(), fontsize=20)
    axs.set_ylabel("PDF", fontsize=20)
    axs.tick_params(width=1.0, labelsize=16)
    if directory is None:
        fig.savefig(f"{key}.png")
    else:
        fig.savefig(f"{directory}/{key}.png")


def plot(
    jd: List[float], rmsd: List[float], title: str, outdir: str, filename: str
) -> None:
    """Plots the Jaccard distance and the RMSD of the pocket heavy atoms. The
    plot will be saved to the specified directory.

    Args:
        jd (list): list containing the Jaccard distance.
        rmsd (list): list containing the RMSD of the pocket heavy atoms.
        title (str): title of the plot.
        outdir (str): path to outdir.
        filename (str): name of the file.
    """

    fig, axs = plt.subplots(1, 2, figsize=(16, 6))
    fig.suptitle(title, fontsize=32)
    axs[0].set_xlabel("Iteration", fontsize=28)
    axs[1].set_xlabel("Iteration", fontsize=28)
    axs[0].set_ylabel("Jaccard distance", fontsize=28)
    axs[1].set_ylabel("Pocket HA RMSD", fontsize=28)
    axs[0].plot(np.arange(0.5, 0.5 * len(jd) + 0.5, 0.5), jd, label="Jaccard", lw=3.0)
    axs[1].plot(
        np.arange(0.5, 0.5 * len(jd) + 0.5, 0.5),
        rmsd,
        label="RMSD",
        lw=3.0,
        color="orange",
    )
    axs[0].tick_params(width=1.0, labelsize=16)
    axs[1].tick_params(width=1.0, labelsize=16)
    plt.savefig(f"{outdir}/{filename}")


def get_pcoords_for_plotting(westfile: str, settings: dict) -> dict:
    """Gets the pcoord and auxdata from the west.h5 file.

    Args:
        westfile (str): path to the west.h5 file.
        settings (dict): dictionary containing the settings for the analysis.

    Returns:
        results (dict): dictionary containing the pcoord and auxdata.
    """

    results = {i: [] for i in settings["pcoord"]}
    for i in settings["auxdata"]:
        results[i] = []
    if "composite" in results:
        if "prmsd" not in results:
            results["prmsd"] = []
        if "bb_rmsd" not in results:
            results["bb_rmsd"] = []

    west = h5py.File(westfile, "r")

    for iteration in list(west["iterations"].keys())[:-1]:
        if len(settings["pcoord"]) == 1:
            results[settings["pcoord"]] = results[settings["pcoord"]] + list(
                west["iterations"][iteration]["pcoord"][:, 1:].flatten()
            )
            for j in settings["auxdata"]:
                results[j] = results[j] + list(
                    west["iterations"][iteration]["auxdata"][j][:, 1:].flatten()
                )

    return results


def plot_stuff(
    jd: List[float],
    min_md_jd: float,
    max_md_jd: float,
    average_jd: float,
    jd_max: List[float],
    jd_min: List[float],
    jd_sd: List[float],
) -> None:
    """Plots the Jaccard distance against iteration. The plot will be saved to
    the specified directory.

    Args:
        jd (list): list containing the Jaccard distance.
        min_md_jd (float): minimum Jaccard distance of the MD simulations.
        max_md_jd (float): maximum Jaccard distance of the MD simulations.
        average_jd (float): average Jaccard distance of the MD simulations.
        jd_max (list): list containing the maximum Jaccard distance.
        jd_min (list): list containing the minimum Jaccard distance.
        jd_sd (list): list containing the standard deviation of the Jaccard distance.
    """

    # TODO: Erich - this plot needs to  be redone
    fig = plt.subplots(1, 1, figsize=(10, 8))
    plt.title("Jaccard distance against iteration", fontsize=20)
    plt.ylabel("JD", fontsize=18)
    plt.errorbar(
        range(len(jd)),
        jd,
        yerr=jd_sd,
        ls="--",
        marker="o",
        capsize=3,
        capthick=1,
        ecolor="gray",
        color="black",
        label="Jaccard distance",
    )
    plt.plot(
        range(len(jd)),
        rolling_average(jd, 5),
        color="red",
        alpha=0.7,
        linewidth=3,
        label="JD rolling average",
    )
    plt.plot(
        range(len(jd)),
        rolling_average(jd_max, 5),
        color="blue",
        alpha=0.7,
        linewidth=3,
        label="Max value rolling average",
    )
    plt.plot(
        range(len(jd)),
        rolling_average(jd_min, 5),
        color="purple",
        alpha=0.7,
        linewidth=3,
        label="Min value rolling average",
    )
    plt.hlines(min_md_jd, 0, len(jd), color="purple", linestyles="dashed", alpha=0.9)
    plt.hlines(max_md_jd, 0, len(jd), color="blue", linestyles="dashed", alpha=0.9)
    plt.hlines(
        average_jd,
        0,
        len(jd),
        color="black",
        linestyles="dashed",
        alpha=0.9,
        label="MD simulation average",
    )
    plt.legend(loc="upper left", bbox_to_anchor=(0, 0.95, 0, 0))


def obtain_paths_highest_rmsd_jd_weight(
    max_rmsd: List[float],
    max_jd: List[float],
    max_weight: List[float],
    reverse_iterations: List[int],
) -> dict:
    """Determines the paths with highest RMSD, JD, and weight.

    Args:
        max_rmsd (list): list containing the maximum RMSD.
        max_jd (list): list containing the maximum Jaccard distance.
        max_weight (list): list containing the maximum weight.
        reverse_iterations (list): list containing the iterations in reverse
            order. (?)

    Returns:
        dict: dictionary containing the highest rmsd, jd, and weight paths.
    """

    # TODO: Erich - I may need to erase this plot

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
        rmsd["jd"].append(
            west["iterations"][iteration]["pcoord"][rmsd["trajectory"][-1]][-1][0]
        )
        rmsd["rmsd"].append(
            west["iterations"][iteration]["pcoord"][rmsd["trajectory"][-1]][-1][1]
        )
        rmsd["trajectory"].append(
            west["iterations"][iteration]["seg_index"][rmsd["trajectory"][-1]][1]
        )
        # obtain progress coordinate and parent walker for walker with highest jaccard distance
        jd["jd"].append(
            west["iterations"][iteration]["pcoord"][jd["trajectory"][-1]][-1][0]
        )
        jd["rmsd"].append(
            west["iterations"][iteration]["pcoord"][jd["trajectory"][-1]][-1][1]
        )
        jd["trajectory"].append(
            west["iterations"][iteration]["seg_index"][jd["trajectory"][-1]][1]
        )
        # obtain progress coordinate and parent walker for walker with highest weight
        weight["jd"].append(
            west["iterations"][iteration]["pcoord"][weight["trajectory"][-1]][-1][0]
        )
        weight["rmsd"].append(
            west["iterations"][iteration]["pcoord"][weight["trajectory"][-1]][-1][1]
        )
        weight["trajectory"].append(
            west["iterations"][iteration]["seg_index"][weight["trajectory"][-1]][1]
        )

    results = {}
    results["max_rmsd"] = rmsd
    results["max_jaccard"] = jd
    results["max_weight"] = weight
    for i in results:
        results[i]["jd"].reverse()
        results[i]["rmsd"].reverse()

    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Script that will read the h5 file and create plots to visualize the progress of the simulation."
    )
    parser.add_argument(
        "h5file",
        type=str,
        help="Define the h5 file file with the settings. It is required",
    )
    parser.add_argument(
        "settings",
        type=str,
        help="Define the json file with the settings. It is required",
    )
    parser.add_argument(
        "outdir",
        type=str,
        help="Define the output directory to put the files in, it will creat it if it does not exist. It is required",
    )
    args = parser.parse_args()

    settings = check_input_settings_file(args.settings)

    # read h5 file and create outdir
    west = h5py.File(args.h5file, "r")

    # make sure outdir exists if not, creat it
    if not os.path.isdir(args.outdir):
        os.system(f"mkdir {args.outdir}")

    # plot num of walkers vs generation
    if settings["subpex"]["plotting"]["num_walkers"]:
        plot_num_walkers(west, args.outdir)
