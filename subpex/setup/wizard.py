# The programs helps set up the SubPEx run. Runs only on linux.

from typing import Tuple

import glob
import os
import re
import sys

import yaml
from loguru import logger

CUR_PATH = os.path.abspath(os.path.dirname(__file__))
PYTHON_EXE = os.path.abspath(sys.executable)
YAML_LOADER = yaml.SafeLoader
STEP_NUM = 1
TEST_MODE = False


def get_step_num() -> str:
    """Returns current step number.

    Returns:
        str: Current step number.
    """

    global STEP_NUM
    STEP_NUM = STEP_NUM + 1
    return str(STEP_NUM - 1)


def get_prelim_sim_files(engine: str) -> None:
    """Asks user for preliminary (equilibrated) simulation files.

    Args:
        engine (str): MD engine to use. "amber" or "namd"
    """

    logger.info(f"STEP {get_step_num()}: Preliminary MD files", "-")
    logger.info(
        "SubPEx assumes you have already run preliminary simulations to equilibrate your system. Where are the coordinate and restart files associated with these preliminary simulations?"
    )

    logger.info(
        "(Note that you can press Ctrl-Z to pause the wizard and search for the file. Type `fg` to resume the wizard when done.)"
    )

    if TEST_MODE:
        logger.info(
            "Downloading Amber .nc, .rst, and .prmtop files for testing (test mode)."
        )
        coor_file = download_testing_files("mol.nc")
        restart_files = [download_testing_files("mol.rst")]
        prmtop_file = download_testing_files("mol.prmtop")
        enter_to_continue()
    else:
        if engine == "amber":
            coor_file = get_choice(
                "coor_file", lambda: file_path("Amber .nc file", ".nc")
            )
            restart_files = [
                get_choice(
                    "restart_files", lambda: file_path("Amber .rst file", ".rst")
                )
            ]
        elif engine == "namd":
            coor_file = get_choice(
                "coor_file", lambda: file_path("NAMD .dcd file", ".dcd")
            )

            restart_files = [
                get_choice("_ignore", lambda: file_path(f"NAMD {ext} file", ext))
                for ext in [".xsc", ".coor", ".inpcrd"]
            ]

            # save_choice("restart_files", restart_files)

        prmtop_file = get_choice(
            "prmtop_file", lambda: file_path("Parameter (.prmtop) file", ".prmtop")
        )

    run_cmd(["rm", "-f", f"{CUR_PATH}/reference/mol.*"])
    link_to_reference_mol(coor_file)
    for restart_file in restart_files:
        link_to_reference_mol(restart_file)
    link_to_reference_mol(prmtop_file)
    clear()


def get_sim_last_frame() -> None:
    """Asks user for the last frame of the preliminary simulation."""

    logger.info(f"STEP {get_step_num()}: Preliminary MD last frame", "-")
    logger.info(
        "SubPEx needs the last frame of the preliminary trajectory as a `pdb` file. You can use programs such as VMD to save the last frame."
    )

    if TEST_MODE:
        logger.info("Downloading Amber .pdb file for testing (test mode).")
        last_frame_file = download_testing_files("mol.pdb")
        enter_to_continue()
    else:
        last_frame_file = get_choice(
            "last_frame_file", lambda: file_path("Last-frame .pdb file", ".pdb")
        )

    link_to_reference_mol(last_frame_file)
    clear()


def update_westcfg_path_vars(westcfg: str) -> str:
    """Updates the path variables in the west.cfg file.

    Args:
        westcfg (str): contents of the west.cfg file as a string.

    Returns:
        str: updated west.cfg file as a string.
    """

    westcfg = re.sub(
        r"\bselection_file: .*?$",
        f"selection_file: {CUR_PATH}/reference/selection_string.txt",
        westcfg,
        flags=re.MULTILINE,
    )

    westcfg = re.sub(
        r"\btopology: .*$",
        f"topology: {CUR_PATH}/reference/mol.prmtop",
        westcfg,
        flags=re.MULTILINE,
    )

    westcfg = re.sub(
        r"\bwest_home: .*$", f"west_home: {CUR_PATH}", westcfg, flags=re.MULTILINE
    )

    westcfg = re.sub(
        r"\breference: .*$",
        f"reference: {CUR_PATH}/reference/mol.pdb",
        westcfg,
        flags=re.MULTILINE,
    )

    westcfg = re.sub(
        r"\breference_fop: .*$",
        f"reference_fop: {CUR_PATH}/reference/fop_ref.xyz",
        westcfg,
        flags=re.MULTILINE,
    )

    return westcfg


def pick_pcoord(westcfg: str) -> Tuple[str, str]:
    """Asks user to pick a progress coordinate.

    Args:
        westcfg (str): contents of the west.cfg file as a string.

    Returns:
        Tuple[str, str]: Tuple containing the updated west.cfg file as a string
        and the progress coordinate name.
    """

    logger.info(f"STEP {get_step_num()}: Progress coordinate", "-")
    logger.info("Which progress coordinate would you like to use?")
    logger.info(
        "composite: combination of pocket and back-bone RMSD (recommended)\nprmsd: RMSD of pocket-lining atoms (not officially supported)\nbb: RMSD of all back-bone atoms (not officially supported)\njd: Jaccard distance between pocket shapes (not officially supported)"
    )

    if TEST_MODE:
        logger.info("Using composite progress coordinate (test mode).")
        pcoord = "composite"
        enter_to_continue()
    else:
        pcoord = get_choice(
            "pcoord",
            lambda: choice(
                "Which progress coordinate?",
                choices=["composite", "prmsd", "bb", "jd", ""],
            ),
        )

    westcfg = re.sub(
        r"(\bpcoord:\s*\n\s*- ).*$", r"\1" + pcoord, westcfg, flags=re.MULTILINE
    )

    clear()
    return westcfg, pcoord


def specify_number_iterations(westcfg: str) -> str:
    """Asks user to specify the number of iterations.

    Args:
        westcfg (str): contents of the west.cfg file as a string.

    Returns:
        str: updated west.cfg file as a string.
    """

    logger.info(f"STEP {get_step_num()}: Number of iterations (generations)", "-")
    logger.info(
        "SubPEx periodically prunes and duplicates short simulations (walkers) to encourage sampling. How many times (iterations/generations) should SubPEx perform this pruning/duplication before stopping? (Recommendation: 30)"
    )

    if TEST_MODE:
        logger.info("Using 3 iterations (test mode).")
        iterations = 3
        enter_to_continue()
    else:
        iterations = int(
            get_choice("iterations", lambda: get_number("Number of iterations"))
        )

    westcfg = re.sub(
        r"\bmax_total_iterations: .*$",
        f"max_total_iterations: {iterations}",
        westcfg,
        flags=re.MULTILINE,
    )

    clear()
    return westcfg


def define_pocket(westcfg: str) -> str:
    """Asks user to define the binding pocket.

    Args:
        westcfg (str): contents of the west.cfg file as a string.

    Returns:
        str: updated west.cfg file as a string.
    """

    logger.info(f"STEP {get_step_num()}: Define the binding pocket", "-")
    logger.info(
        "You must define the location of the binding pocket you wish to sample. Visual inspection (e.g., using VMD) is often useful for this step. See README.md for more suggestions."
    )

    logger.info(
        f"Be sure to specify the location of the pocket relative to the reference PDB file: {CUR_PATH}/reference/mol.pdb"
    )

    if TEST_MODE:
        logger.info(
            "Using default pocket: radius = 6.5, center = [30.0, 41.5, 30.4] (test mode)."
        )
        pocket_radius = 6.5
        x_coor = 30.0
        y_coor = 41.5
        z_coor = 30.4
        enter_to_continue()
    else:
        pocket_radius = get_choice("pocket_radius", lambda: get_number("Pocket radius"))
        x_coor = get_choice(
            "x_coor", lambda: get_number("X coordinate of pocket center")
        )
        y_coor = get_choice(
            "y_coor", lambda: get_number("Y coordinate of pocket center")
        )
        z_coor = get_choice(
            "z_coor", lambda: get_number("Z coordinate of pocket center")
        )

    westcfg = re.sub(
        r"\bcenter: \[.*\]",
        f"center: [{str(x_coor)}, {str(y_coor)}, {str(z_coor)}]",
        westcfg,
        flags=re.MULTILINE,
    )

    westcfg = re.sub(
        r"\bradius: .*$", f"radius: {str(pocket_radius)}", westcfg, flags=re.MULTILINE
    )

    with open("west.cfg", "w") as f:
        f.write(westcfg)

    run_cmd(["rm", "-f", f"{CUR_PATH}/reference/fop_ref.xyz"])
    run_cmd(["rm", "-f", f"{CUR_PATH}/reference/selection_string.txt"])
    run_cmd([PYTHON_EXE, "./westpa_scripts/get_reference_fop.py", "west.cfg"])
    if not os.path.exists(f"{CUR_PATH}/reference/fop_ref.xyz") or not os.path.exists(
        f"{CUR_PATH}/reference/selection_string.txt"
    ):
        logger.info(
            (
                (
                    "ERROR: Could not generate reference-fop file (fop_ref.xyz) and/or selection-string file (selection_string.txt) in the directory "
                    + CUR_PATH
                )
                + "/reference/"
            )
        )

        logger.info(
            "Are you sure the center and radius specified encompass at least some of the pocket-lining residues?"
        )

        sys.exit(1)
    clear()
    westcfg = check_pocket(westcfg)
    return westcfg


def check_pocket(westcfg: str) -> str:
    """Asks user to check the binding pocket.

    Args:
        westcfg (str): contents of the west.cfg file as a string.

    Returns:
        str: updated west.cfg file as a string.
    """

    logger.info(f"STEP {get_step_num()}: Check the binding pocket", "-")
    logger.info(
        f"The wizard has generated a pocket-field-of-points (fop) file ({CUR_PATH}/reference/fop_ref.xyz) and a selection-string file ({CUR_PATH}/reference/selection_string.txt)."
    )

    logger.info("Please visually inspect these files to:")
    logger.info(
        "1. Ensure that the points in fop_ref.xyz entirely fill the pocket of interest."
    )

    logger.info(
        "2. Ensure that the residues described in selection_string.txt truly line the pocket of interest."
    )

    logger.info(
        "(Note that the popular molecular visualization program VMD can load xyz files and select residues."
    )

    if TEST_MODE:
        logger.info("Assuming pocket correct (test mode).")
        enter_to_continue()
    elif choice("Can you confirm the pocket is properly defined?") == "n":
        clear()
        forget_choice("pocket_radius")
        forget_choice("x_coor")
        forget_choice("y_coor")
        forget_choice("z_coor")
        westcfg = define_pocket(westcfg)
    clear()
    return westcfg


def update_envsh(envsh: str, engine: str) -> str:
    """Asks user to update the env.sh file.

    Args:
        envsh (str): contents of the env.sh file as a string.
        engine (str): engine to use.

    Returns:
        str: updated env.sh file as a string.
    """

    logger.info(f"STEP {get_step_num()}: Update env.sh", "-")
    envsh = re.sub(
        r'^export ENGINE=".*?"',
        f'export ENGINE="{engine.upper()}"',
        envsh,
        flags=re.MULTILINE,
    )

    if TEST_MODE:
        logger.info("Using GPU (test mode).")
        mode = "GPU"
        enter_to_continue()
    else:
        mode = get_choice(
            "mode",
            lambda: choice(
                "How will you run SubPEx (only GPU officially supported)?",
                choices=["GPU", "MPI", "MULTITHREAD"],
            ),
        )

    envsh = re.sub(
        r'^export MODE=".*?"',
        f'export MODE="{mode.upper()}"',
        envsh,
        flags=re.MULTILINE,
    )

    clear()
    return envsh


def define_adaptive_bins(adaptivepy: str, pcoord: str) -> str:
    """Asks user to define the adaptive bins.

    Args:
        adaptivepy (str): contents of the adaptive.py file as a string.
        pcoord (str): progress coordinate to use.

    Returns:
        str: updated adaptive.py file as a string.
    """

    # If pcoord is "jd", max is 1.0.
    logger.info(f"STEP {get_step_num()}: Define adaptive bins", "-")
    logger.info(
        "SubPEx uses WESTPA to enhance the sampling of your pocket. A progress coordinate assesses the extent of sampling. This progress coordinate is divided into bins, and WESTPA works to make sure the simulations collectively sample all bins equally."
    )

    logger.info(f"You previously selected the progress coordinate: {pcoord}")
    if pcoord == "jd":
        maxcap = 1.0
        save_choice("maxcap", maxcap)
    else:
        logger.info(
            "What is the maximum value of the progress coordinate beyond which SubPEx should no longer enhance sampling (in Angstroms)? We recommend 5."
        )

        if TEST_MODE:
            logger.info("Using 5.0 (test mode).")
            maxcap = 5.0
            enter_to_continue()
        else:
            maxcap = get_choice("maxcap", lambda: get_number("Maximum value"))

    logger.info("\nHow many bins would you like to use? We recommend 15.")

    if TEST_MODE:
        logger.info("Using 15 bins (test mode).")
        binsperdim = 15
        enter_to_continue()
    else:
        binsperdim = get_choice("binsperdim", lambda: int(get_number("Bin count")))

    logger.info(
        "\nHow many walkers (mini simulations) would you like run per bin? Using more improves sampling at the cost of computer resources. We recommend 3."
    )

    if TEST_MODE:
        logger.info("Using 3 walkers per bin (test mode).")
        bintargetcount = 3
        enter_to_continue()
    else:
        bintargetcount = get_choice(
            "bintargetcount", lambda: int(get_number("Number of walkers per bin"))
        )

    adaptivepy = re.sub(
        r"^maxcap\s*=\s*\[.*\]",
        f"maxcap = [{str(maxcap)}]",
        adaptivepy,
        flags=re.MULTILINE,
    )

    adaptivepy = re.sub(
        r"^binsperdim\s*=\s*\[.*\]",
        f"binsperdim = [{str(binsperdim)}]",
        adaptivepy,
        flags=re.MULTILINE,
    )

    adaptivepy = re.sub(
        r"^bintargetcount\s*=\s*\d*",
        f"bintargetcount = {str(bintargetcount)}",
        adaptivepy,
        flags=re.MULTILINE,
    )

    clear()
    return adaptivepy


def update_run_time_and_job_name(westcfg: str) -> str:
    """Asks user to update the run time and job name.

    Args:
        westcfg (str): contents of the west.cfg file as a string.

    Returns:
        str: updated west.cfg file as a string.
    """

    logger.info(f"STEP {get_step_num()}: Update job run time and name", "-")

    logger.info(
        "How long should your SubPEx job run? (Format your response like HH:MM:SS, e.g., 72:00:00 for 72 hours)"
    )
    if TEST_MODE:
        logger.info("Using 72:00:00 (test mode).")
        run_time = "72:00:00"
        enter_to_continue()
    else:
        run_time = get_choice("run_time", lambda: get_time("Run time"))

    logger.info('\nWhat is the name of your SubPEx job? (e.g., "my_job")')

    if TEST_MODE:
        logger.info('Using "subpex_job" (test mode).')
        job_name = "subpex_job"
        enter_to_continue()
    else:
        job_name = get_choice("job_name", lambda: input("Job name: ").replace('"', ""))

    westcfg = re.sub(
        r"\bmax_run_wallclock:\s*\d{1,2}:\d{1,2}:\d{1,2}\b",
        f"max_run_wallclock:    {run_time}",
        westcfg,
        flags=re.MULTILINE,
    )

    for submit_file in glob.glob("aux_scripts/run.slurm*.sh"):
        with open(submit_file) as f:
            runsh = f.read()
        runsh = re.sub(
            r"^#SBATCH --time=.*$",
            f"#SBATCH --time={run_time}",
            runsh,
            flags=re.MULTILINE,
        )
        runsh = re.sub(
            r"^#SBATCH --job-name=.*$",
            f"#SBATCH --job-name={job_name}",
            runsh,
            flags=re.MULTILINE,
        )
        with open(submit_file, "w", encoding="utf-8") as f:
            f.write(runsh)

    clear()
    return westcfg


def run_init() -> str:
    """Runs the init.sh script to initialize the SubPEx run.

    Returns:
        str: contents of the west.cfg file as a string.
    """

    logger.info(f"STEP {get_step_num()}: Initialize the SubPEx run", "-")
    if os.path.exists("west.h5"):
        logger.info(
            "Would you like to clear the previous SubPEx run? WARNING: This will delete previously generated files, erasing the previous results so you can start a fresh SubPEx run."
        )

        if choice("Delete previous run and start over") == "y":
            print("")
            run_cmd(["rm", "-f", "./job_logs/*"])
            run_cmd(["rm", "-f", "west.h5"])
            run_cmd(["./init.sh"])
            if not os.path.exists("west.h5"):
                logger.info("ERROR: Could not initialize the SubPEx run.")
                sys.exit(1)
            else:
                clear()
    else:
        run_cmd(["rm", "-f", "./job_logs/*"])
        run_cmd(["./init.sh"])
        if not os.path.exists("west.h5"):
            logger.info("ERROR: Could not initialize the SubPEx run.")
            sys.exit(1)
        else:
            clear()
    return westcfg


def run_wizard():
    # Load files
    with open("./west.cfg", "r", encoding="utf-8") as f:
        westcfg = f.read()

    with open("./adaptive_binning/adaptive.py", "r", encoding="utf-8") as f:
        adaptivepy = f.read()

    with open("./westpa_scripts/get_pcoord.sh", "r", encoding="utf-8") as f:
        get_pcoord = f.read()

    with open("./westpa_scripts/runseg.sh", "r", encoding="utf-8") as f:
        runseg = f.read()

    with open("./env.sh", "r", encoding="utf-8") as f:
        envsh = f.read()

    if check_if_restart_sim():
        # Exit
        sys.exit(0)

    check_existing_params()
    confirm_environment()
    confirm_dependencies()
    engine, get_pcoord, runseg = amber_or_namd(get_pcoord, runseg)
    get_prelim_sim_files(engine)
    get_sim_last_frame()
    westcfg = update_westcfg_path_vars(westcfg)
    westcfg, pcoord = pick_pcoord(westcfg)
    westcfg = specify_number_iterations(westcfg)
    westcfg = define_pocket(westcfg)
    envsh = update_envsh(envsh, engine)
    adaptivepy = define_adaptive_bins(adaptivepy, pcoord)
    westcfg = update_run_time_and_job_name(westcfg)

    # Save west.cfg
    with open("./west.cfg", "w", encoding="utf-8") as f:
        f.write(westcfg)

    with open("./adaptive_binning/adaptive.py", "w", encoding="utf-8") as f:
        f.write(adaptivepy)

    with open("./westpa_scripts/get_pcoord.sh", "w", encoding="utf-8") as f:
        f.write(get_pcoord)

    with open("./westpa_scripts/runseg.sh", "w", encoding="utf-8") as f:
        f.write(runseg)

    with open("./env.sh", "w", encoding="utf-8") as f:
        f.write(envsh)

    run_init()
    finished()
