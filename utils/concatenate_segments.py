import MDAnalysis as mda
import glob
import os
import argparse
import json


def concatenate_segments(directory, prmtop, reference, filename="concatenated", selection=None, outdir="concat", stride=False, num_frames=None):
    """

    :param directory:
    :param prmtop:
    :param reference:
    :param filename:
    :param selection:
    :param outdir:
    :param stride:
    :param num_frames:
    :return:
    """
    os.system("mkdir {}".format(outdir))
    all_iters = []
    for iteration in sorted(glob.glob(directory + "traj_segs/*/")):
        walkers = sorted(glob.glob(iteration + "/*/"))
        dcds = ""
        for i in walkers:
            dcds += " {}seg.dcd".format(i)
        if stride:
            command = "prody catdcd -o {}/{}/{}.dcd --psf {} --stride {} --first {}".format(
                dir=directory,
                outdir=outdir,
                filename=iteration.split("/")[-2],
                prmtop=prmtop,
                num_frames=num_frames,
                last_frame=num_frames - 1)

        else:
            command = "prody catdcd -o {dir}/{outdir}/{filename}.dcd ".format(
                dir=directory,
                outdir=outdir,
                filename=iteration.split("/")[-2],
                prmtop=prmtop)
        command += dcds
        os.system(command)
        all_iters.append("{dir}/{outdir}/{filename}.dcd".format(dir=directory,
                                                                outdir=outdir,
                                                                filename=iteration.split("/")[-2]))

    command = "prody catdcd -o {dir}/{outdir}/{filename}.dcd ".format(dir=directory,
                                                                      outdir=outdir,
                                                                      filename=filename)
    for i in all_iters:
        command += " "
        command += i
    os.system(command)
    for i in all_iters:
        os.system("rm " + i)

    universe = mda.Universe(prmtop, "{dir}/{outdir}/{filename}.dcd".format(dir=directory,
                                                                           outdir=outdir,
                                                                           filename=filename))
    with open(selection, "r") as f:
        selection_string = f.readlines()[0]

    pocket = universe.select_atoms(selection_string)
    pocket.write("{dir}/{outdir}/{filename}_pocket_only.dcd".format(dir=directory,
                                                                    outdir=outdir,
                                                                    filename=filename), frames="all")
    pocket.write("{dir}/{outdir}/{filename}_pocket_only.pdb".format(dir=directory,
                                                        outdir=outdir,
                                                        filename=filename))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Script that will use prody catdcd and mdanalysis to concatenate seg.dcd files.")
    parser.add_argument("settings", type=str, help="Define the json file with the settings. It is required")
    parser.add_argument("home_dir", type=str, help="Define the home directory of the SubPEx run. It is required")
    parser.add_argument("filename", type=str, help="Define the filename for the concatenated files. It is required")
    parser.add_argument("--stride", action="store_true")
    parser.add_argument("--outdir", type=str, help="will save results in specified directory")
    args = parser.parse_args()

    try:
        with open(args.settings, "r") as f:
            settings = json.load(f)
    except IOError:
        print("Could not load the json file with the settings")
        print("make sure the file exists and is correctly formatted")
        raise IOError("Could not load the json file with the settings")

    if args.stride:
        if settings["calculated_points"] == -1:
            segments = glob.glob("{dir}/{traj_segs/*/*/seg.dcd")
            universe = mda.Universe(settings["topology"], segments[0])
            num_frames = len(universe.trajectory)
        else:
            num_frames = settings["calculated_points"]
    else:
        num_frames = None
    if args.outdir is None:
        outdir = "concatenated"
    else:
        outdir = args.outdir
    concatenate_segments(args.home_dir + "/", settings["topology"], settings["reference"], filename=args.filename,
                         selection=settings["selection_file"], outdir=outdir, stride=False, num_frames=None)
