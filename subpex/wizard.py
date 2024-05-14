# The programs helps set up the SubPEx run. Runs only on linux.

from typing import Any, List, Tuple, Union

import glob
import json
import os
import re
import subprocess
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


def choice(prompt: str, choices: List[str] | None = None) -> str:
    """Prompts user for input, with optional choices.

    Args:
        prompt: Prompt to display.
        choices: List of choices. Defaults to None.

    Returns:
        str: User's choice.
    """

    if choices is None:
        choices = ["y", "n"]

    lowers = [c.lower() for c in choices]
    choices_to_show = [c for c in choices if c != ""]
    while True:
        answer = input(f"{prompt} (" + "/".join(choices_to_show) + "): ").lower()
        if answer not in lowers:
            if len(choices_to_show) == 2:
                print("Please answer " + " or ".join(choices_to_show))
            else:
                print(
                    "Please answer "
                    + ", ".join(choices_to_show[:-1])
                    + ", or "
                    + choices_to_show[-1]
                )
            continue
        return answer


def clear():
    """Clears screen."""
    os.system("clear")


def file_path(prompt: str, ext: str) -> str:
    """Prompts user for file path, with optional extension.

    Args:
        prompt (str): Prompt to display.
        ext (str): File extension (for checking).

    Returns:
        str: User's choice.
    """

    while True:
        answer = input(f"{prompt}: ").strip()
        pth = os.path.abspath(answer)
        if not os.path.exists(pth):
            logger.info("\nFile does not exist: " + pth + ". Please try again.")
            continue
        if not pth.lower().endswith(ext.lower()):
            logger.info("\nFile does not end in ." + ext + ". Please try again.")
            continue
        return pth


def get_time(prompt: str) -> str:
    """Prompts user for time in format HH:MM:SS.

    Args:
        prompt (str): Prompt to display.

    Returns:
        str: User's choice.
    """

    while True:
        answer = input(f"{prompt}: ").strip()

        # Make sure answer matches format HH:MM:SS
        if not re.match(r"^\d{1,2}:\d{1,2}:\d{1,2}$", answer):
            logger.info("\nPlease enter time in the format HH:MM:SS. Try again.")
            continue

        return answer


def save_choice(key: str, val: Any) -> None:
    """Saves choice to wizard.saved file.

    Args:
        key (str): Key to save.
        val (any): Value to save.
    """

    if os.path.exists("wizard.saved"):
        with open("wizard.saved", "r") as f:
            config = json.load(f)
    else:
        config = {}
    config[key] = val
    with open("wizard.saved", "w") as f:
        json.dump(config, f, indent=4)


def get_choice(key: str, func_if_absent: callable) -> Any:
    """Gets choice from wizard.saved file, or prompts user if not present.

    Args:
        key (str): Key to get.
        func_if_absent (callable): Function to call if key is not present.

    Returns:
        any: User's choice.
    """

    if not os.path.exists("wizard.saved"):
        val = func_if_absent()
        save_choice(key, val)
        return val
    with open("wizard.saved", "r") as f:
        config = json.load(f)
    if key not in config:
        val = func_if_absent()
        save_choice(key, val)
        return val
    logger.info(f"Using previous choice for {key}: {str(config[key])}")
    return config[key]


def forget_choice(key: str):
    """Forgets choice from wizard.saved file.

    Args:
        key (str): Key to forget.
    """

    if os.path.exists("wizard.saved"):
        with open("wizard.saved", "r") as f:
            config = json.load(f)
    else:
        config = {}

    if key in config:
        del config[key]

    with open("wizard.saved", "w") as f:
        json.dump(config, f, indent=4)


def run_cmd(prts: List[str]) -> str:
    """Runs command, prints output, and returns output.

    Args:
        prts: Command to run.

    Returns:
        Output of command.
    """

    # if prts is a list
    if isinstance(prts, list):
        logger.info("Running command: " + " ".join(prts))
        shell = False
    else:
        logger.info(f"Running command: {prts}")
        shell = True

    def run_iter(prts):
        cmd = subprocess.Popen(
            prts,
            stdout=subprocess.PIPE,
            universal_newlines=True,
            bufsize=1,
            shell=shell,
        )
        yield from iter(cmd.stdout.readline, "")
        cmd.stdout.close()
        return_code = cmd.wait()
        if return_code:
            error = cmd.stderr.read()
            print(error)
            raise subprocess.CalledProcessError(return_code, cmd)

    output = ""
    for line in run_iter(prts):
        print(line, end="")
        output += line

    return output.strip()


def check_if_restart_sim() -> bool:
    """Checks if user wants to restart simulation.

    Returns:
        bool: True if user wants to restart simulation, False otherwise.
    """

    if os.path.exists("west.h5"):
        logger.info(f"STEP {get_step_num()}: Restart previous SubPEx run?", "-")
        logger.info(
            "It appears you have previously run SubPEx in this directory because a west.h5 exists."
        )
        logger.info(
            "Would you like to have SubPEx resume the previous run rather than start over?"
        )
        if choice("Resume previous run?") == "y":
            run_cmd("sp_restart -n $(ls traj_segs/ | sort -n | tail -n 1)")
            logger.info(
                "\nYour SubPEx job is now ready for a restart run, using ./run.sh."
            )
            # clear()
            return True

    clear()
    return False


def check_existing_params():
    """Checks if user wants to use existing parameters."""

    if not os.path.exists("wizard.saved"):
        return

    logger.info(f"STEP {get_step_num()}: Use previous wizard choices?", "-")
    logger.info(
        'You previously used this wizard, and your choices were saved to the file "wizard.saved". Here are the previous choices:'
    )

    with open("wizard.saved", "r") as f:
        config = json.load(f)
    for key in config:
        print(f"  {key}: {config[key]}")

    logger.info("\nWould you like to use the same choices this time?")

    if choice("Use previous choices?") == "n":
        os.system("rm wizard.saved")
    clear()


def confirm_environment():
    """Confirms that user has created and activated anaconda environment."""

    logger.info(f"STEP {get_step_num()}: Environment", "-")
    logger.info(
        "You should run this wizard using an appropriate anaconda environment. To create the environment, install anaconda and create/activate the environment like this:"
    )

    logger.info("conda env create -f environment.yaml\n conda activate westpa")
    if choice("Have you created and activated the environment?") == "n":
        logger.info("\nPlease create and activate the environment. Then try again.")
        exit(1)
    clear()


def confirm_dependencies():
    """Confirms that user has installed dependencies."""

    logger.info(f"Step {get_step_num()}: Check dependencies", "-")

    if (
        choice("Would you like to check that key dependencies have been installed?")
        == "y"
    ):
        print("")

        # Load environment.yaml into a dictionary
        with open("environment.yaml", "r") as f:
            env = yaml.load(f, Loader=YAML_LOADER)
        env = [v for v in env["dependencies"] if v[:1] != "_"]
        env = [v.split("=")[0] for v in env]
        to_ignore = [
            "biopython",
            "bzip2",
            "curl",
            "freetype",
            "griddataformats",
            "hdf4",
            "hdf5",
            "ipython",
            "jpeg",
            "krb5",
            "lcms2",
            "libblas",
            "libcblas",
            "libcurl",
            "libedit",
            "libev",
            "libffi",
            "libgfortran5",
            "libgomp",
            "liblapack",
            "libnetcdf",
            "libnghttp2",
            "libopenblas",
            "libpng",
            "libsodium",
            "libssh2",
            "libtiff",
            "mdanalysis",
            "ncurses",
            "netcdf4",
            "openssl",
            "pillow",
            "python",
            "python_abi",
            "pyyaml",
            "pyzmq",
            "sqlite",
            "tk",
            "westpa",
            "xz",
            "zeromq",
            "zstd",
        ]

        for v in env:
            if v in to_ignore:
                continue
            if "-" in v:
                continue
            print(f"Testing module {v}")
            try:
                __import__(v)
            except Exception:
                logger.info(
                    (
                        "\nModule "
                        + v
                        + " is missing. Have you created and activated the westpa environment?"
                    )
                )

                logger.info(
                    "conda env create -f environment.yaml\n conda activate westpa"
                )
                logger.info("Please create/activate it and try again.")
                exit(1)
    clear()


def amber_or_namd(get_pcoord: str, runseg: str) -> Tuple[str, str, str]:
    """Asks user if they want to use Amber or NAMD.

    Args:
        get_pcoord (str): contents of the get_pcoord.sh file as a string.
        runseg (str): contents of the runseg.sh file as a string.

    Returns:
        Tuple[str, str, str]: Tuple containing the engine, get_pcoord.sh file
        as a string, and runseg.sh file as a string.
    """

    logger.info(f"STEP {get_step_num()}: MD engine", "-")

    if TEST_MODE:
        engine = "amber"
        logger.info("Using Amber MD engine (test mode).")
        enter_to_continue()
    else:
        engine = get_choice(
            "engine",
            lambda: choice(
                "Will you use the Amber or NAMD MD engine?", choices=["Amber", "NAMD"]
            ),
        )

    clear()
    get_pcoord = re.sub(
        r'^export ENGINE=".*?"',
        f'export ENGINE="{engine.upper()}"',
        get_pcoord,
        flags=re.MULTILINE,
    )

    runseg = re.sub(
        r'^export ENGINE=".*"',
        f'export ENGINE="{engine.upper()}"',
        runseg,
        flags=re.MULTILINE,
    )

    return engine, get_pcoord, runseg


def download_testing_files(filename: str) -> str:
    """Downloads testing files from durrantlab.

    Args:
        filename (str): Name of file to download.

    Returns:
        str: Path to downloaded file.
    """

    url = f"https://durrantlab.pitt.edu/apps/subpex/{filename}"

    # Make the downloads directory if it doesn't exist
    if not os.path.exists("downloads"):
        os.mkdir("downloads")

    # Download the file using a library from the python standard library
    import urllib.request

    logger.info(f"Downloading {filename} from {url}...")

    urllib.request.urlretrieve(url, f"downloads/{filename}")

    return os.path.abspath(f"downloads/{filename}")


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
