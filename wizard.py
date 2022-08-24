# The programs helps set up the SubPEx run. Runs only on linux.

import textwrap
import os
import json
import sys
import subprocess
import time
import re
import glob

try:
    import yaml
except:
    print("\nPyYAML is missing. Have you set up your westpa conda environment?\n")
    print("conda env create -f environment.yaml\nconda activate westpa\n")
    exit(1)

if os.name == "nt":
    print("\nWizard script supported only on Linux. Sorry!\n")
    exit(1)

# Make sure west.cfg in present directory
if not os.path.exists("west.cfg"):
    print(
        "\nPlease run this script from the same directory containing the west.cfg file.\n"
    )
    exit(1)

CUR_PATH = os.path.abspath(os.path.dirname(__file__))
PYTHON_EXE = os.path.abspath(sys.executable)
YAML_LOADER = yaml.SafeLoader
STEP_NUM = 1


def get_step_num():
    global STEP_NUM
    STEP_NUM = STEP_NUM + 1
    return str(STEP_NUM - 1)


def log(txt, underlined=None):
    txt = txt.split("\n")
    for line in txt:
        t = textwrap.fill(line.strip(), 80)
        print(t)
        if underlined:
            print(underlined * len(t))
    print("")


def choice(prompt, choices=None):
    if choices is None:
        choices = ["y", "n"]

    lowers = [c.lower() for c in choices]
    while True:
        answer = input(f"{prompt} (" + "/".join(choices) + "): ").lower()
        if answer not in lowers:
            if len(choices) == 2:
                print("Please answer " + " or ".join(choices))
            else:
                print(
                    "Please answer " + ", ".join(choices[:-1]) + ", or " + choices[-1]
                )
            continue
        return answer


def clear():
    os.system("clear")


def file_path(prompt, ext):
    while True:
        answer = input(f"{prompt}: ").strip()
        pth = os.path.abspath(answer)
        if not os.path.exists(pth):
            log("\nFile does not exist: " + pth + ". Please try again.")
            continue
        if not pth.lower().endswith(ext.lower()):
            log("\nFile does not end in ." + ext + ". Please try again.")
            continue
        return pth


def get_number(prompt):
    while True:
        answer = input(f"{prompt}: ").strip()
        try:
            answer = float(answer)
        except:
            log("\nPlease enter a number. Try again.")
            continue
        return answer


def get_time(prompt):
    while True:
        answer = input(f"{prompt}: ").strip()

        # Make sure answer matches format HH:MM:SS
        if not re.match(r"^\d{1,2}:\d{1,2}:\d{1,2}$", answer):
            log("\nPlease enter time in the format HH:MM:SS. Try again.")
            continue

        return answer


def save_choice(key, val):
    if os.path.exists("wizard.saved"):
        with open("wizard.saved", "r") as f:
            config = json.load(f)
    else:
        config = {}
    config[key] = val
    with open("wizard.saved", "w") as f:
        json.dump(config, f, indent=4)


def get_choice(key, func_if_absent):
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
    log(f"Using previous choice for {key}: {str(config[key])}")
    return config[key]


def forget_choice(key):
    if os.path.exists("wizard.saved"):
        with open("wizard.saved", "r") as f:
            config = json.load(f)
    else:
        config = {}

    if key in config:
        del config[key]

    with open("wizard.saved", "w") as f:
        json.dump(config, f, indent=4)


def run_cmd_old(prts):
    log("Running command: " + " ".join(prts))
    cmd = subprocess.Popen(
        prts,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    cmd.wait()

    output, errors = cmd.communicate()
    output = output.strip()
    errors = errors.strip()

    if output:
        log(output)
    if errors:
        log(errors)
    
    return output + "\n\n" + errors


def run_cmd(prts):
    # if prts is a list
    if isinstance(prts, list):
        log("Running command: " + " ".join(prts))
        shell = False
    else:
        log("Running command: " + prts)
        shell = True

    def run_iter(prts):
        cmd = subprocess.Popen(
            prts, stdout=subprocess.PIPE, universal_newlines=True, bufsize=1,
            shell=shell
        )
        for stdout_line in iter(cmd.stdout.readline, ""):
            yield stdout_line
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


def link_to_reference_mol(flnm):
    ext = flnm.split(".")[-1]
    dest_file = f"{CUR_PATH}/reference/mol.{ext}"
    if os.path.exists(dest_file):
        os.remove(dest_file)
    run_cmd(["ln", "-s", flnm, dest_file])


def openfile(filename):
    if not os.path.exists(f"{filename}.orig"):
        os.system(f"cp {filename} {filename}.orig")
    return open(f"{filename}.orig", "r")


def check_if_restart_sim():
    if os.path.exists("west.h5"):
        log(f"STEP {get_step_num()}: Restart previous SubPEx run?", "-")
        log("It appears you have previously run SubPEx in this directory because a west.h5 exists.")
        log("Would you like to have SubPEx resume the previous run rather than start over?")
        if choice("Resume previous run?") == "y":
            run_cmd("./utils/restart_subpex.sh -n $(ls traj_segs/ | sort -n | tail -n 1)")
            log("\nYour SubPEx job is now ready for a restart run, likely using one of the run*.sh files.")
            # clear()
            return True

    clear()
    return False



def check_existing_params():
    if os.path.exists("wizard.saved"):
        log(f"STEP {get_step_num()}: Use previous wizard choices?", "-")
        log(
            "You previously used this wizard, and your choices were saved. Would you like to use the same choices this time?"
        )
        if choice("Use previous choices?") == "n":
            os.system("rm wizard.saved")
        clear()


def confirm_environment():
    log(f"STEP {get_step_num()}: Environment", "-")
    log(
        "You should run this wizard using an appropriate anaconda environment. To create the environment, install anaconda and create/activate the environment like this:"
    )

    log("conda env create -f environment.yaml\n conda activate westpa")
    if choice("Have you created and activated the environment?") == "n":
        log("\nPlease create and activate the environment. Then try again.")
        exit(1)
    clear()


def confirm_dependencies():
    log(f"Step {get_step_num()}: Check dependencies", "-")

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
        except:
            log(
                (
                    "\nModule "
                    + v
                    + " is missing. Have you created and activated the westpa environment?"
                )
            )

            log("conda env create -f environment.yaml\n conda activate westpa")
            log("Please create/activate it and try again.")
            exit(1)
    clear()


def amber_or_namd(get_pcoord, runseg):
    log(f"STEP {get_step_num()}: MD engine", "-")

    engine = get_choice(
        "engine",
        lambda: choice(
            "Will you use the Amber or NAMD MD engine?", choices=["Amber", "NAMD"]
        ),
    )

    clear()
    get_pcoord = re.sub(
        r'^export ENGINE=".*?"',
        (('export ENGINE="' + engine.upper()) + '"'),
        get_pcoord,
        flags=re.MULTILINE,
    )

    runseg = re.sub(
        r'^export ENGINE=".*"',
        'export ENGINE="' + engine.upper() + '"',
        runseg,
        flags=re.MULTILINE,
    )

    return engine, get_pcoord, runseg


def get_prelim_sim_files(engine):
    log(f"STEP {get_step_num()}: Preliminary MD files", "-")
    log(
        "SubPEx assumes you have already run preliminary simulations to equilibrate your system. Where are the coordinate and restart files associated with these preliminary simulations?"
    )

    log(
        "(Note that you can press Ctrl-Z to pause the wizard and search for the file. Type `fg` to resume the wizard when done.)"
    )

    if engine == "amber":
        coor_file = get_choice("coor_file", lambda: file_path("Amber .nc file", ".nc"))
        restart_files = [
            get_choice("restart_files", lambda: file_path("Amber .rst file", ".rst"))
        ]

    elif engine == "namd":
        coor_file = get_choice("coor_file", lambda: file_path("NAMD .dcd file", ".dcd"))

        restart_files = [
            get_choice("_ignore", lambda: file_path(f"NAMD {ext} file", ext))
            for ext in [".xsc", ".coor", ".inpcrd"]
        ]

        save_choice("restart_files", restart_files)
    prmtop_file = get_choice(
        "prmtop_file", lambda: file_path("Parameter (.prmtop) file", ".prmtop")
    )

    run_cmd(["rm", "-f", f"{CUR_PATH}/reference/mol.*"])
    link_to_reference_mol(coor_file)
    for restart_file in restart_files:
        link_to_reference_mol(restart_file)
    link_to_reference_mol(prmtop_file)
    clear()


def get_sim_last_frame():
    log(f"STEP {get_step_num()}: Preliminary MD last frame", "-")
    log(
        "SubPEx needs the last frame of the preliminary trajectory as a `pdb` file. You can use programs such as VMD to save the last frame."
    )

    last_frame_file = get_choice(
        "last_frame_file", lambda: file_path("Last-frame .pdb file", ".pdb")
    )

    link_to_reference_mol(last_frame_file)
    clear()


def update_westcfg_path_vars(westcfg):
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


def pick_pcoord(westcfg):
    log(f"STEP {get_step_num()}: Progress coordinate", "-")
    log("Which progress coordinate would you like to use?")
    log(
        "composite: combination of pocket and back-bone RMSD (recommended)\nprmsd: RMSD of pocket-lining atoms\nbb: RMSD of all back-bone atoms\njd: Jaccard distance between pocket shapes"
    )

    pcoord = get_choice(
        "pcoord",
        lambda: choice(
            "Which progress coordinate?", choices=["composite", "prmsd", "bb", "jd"]
        ),
    )

    westcfg = re.sub(
        r"(\bpcoord:\s*\n\s*- ).*$", r"\1" + pcoord, westcfg, flags=re.MULTILINE
    )

    clear()
    return westcfg, pcoord


def define_pocket(westcfg):
    log(f"STEP {get_step_num()}: Define the binding pocket", "-")
    log(
        "You must define the location of the binding pocket you wish to sample. Visual inspection (e.g., using VMD) is often useful for this step. See README.md for more suggestions."
    )

    log(
        f"Be sure to specify the location of the pocket relative to the reference PDB file: {CUR_PATH}/reference/mol.pdb"
    )

    pocket_radius = get_choice("pocket_radius", lambda: get_number("Pocket radius"))

    x_coor = get_choice("x_coor", lambda: get_number("X coordinate of pocket center"))

    y_coor = get_choice("y_coor", lambda: get_number("Y coordinate of pocket center"))

    z_coor = get_choice("z_coor", lambda: get_number("Z coordinate of pocket center"))

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
        log(
            (
                (
                    "ERROR: Could not generate reference-fop file (fop_ref.xyz) and/or selection-string file (selection_string.txt) in the directory "
                    + CUR_PATH
                )
                + "/reference/"
            )
        )

        log(
            "Are you sure the center and radius specified encompass at least some of the pocket-lining residues?"
        )

        sys.exit(1)
    clear()
    westcfg = check_pocket(westcfg)
    return westcfg


def check_pocket(westcfg):
    log(f"STEP {get_step_num()}: Check the binding pocket", "-")
    log(
        f"The wizard has generated a pocket-field-of-points (fop) file ({CUR_PATH}/reference/fop_ref.xyz) and a selection-string file ({CUR_PATH}/reference/selection_string.txt)."
    )

    log("Please visually inspect these files to:")
    log(
        "1. Ensure that the points in fop_ref.xyz entirely fill the pocket of interest."
    )

    log(
        "2. Ensure that the residues described in selection_string.txt truly line the pocket of interest."
    )

    log(
        "(Note that the popular molecular visualization program VMD can load xyz files and select residues."
    )

    if choice("Can you confirm the pocket is properly defined?") == "n":
        clear()
        forget_choice("pocket_radius")
        forget_choice("x_coor")
        forget_choice("y_coor")
        forget_choice("z_coor")
        westcfg = define_pocket(westcfg)
    clear()
    return westcfg


def update_envsh(envsh, engine):
    log(f"STEP {get_step_num()}: Update env.sh", "-")
    envsh = re.sub(
        r'^export ENGINE=".*?"',
        (('export ENGINE="' + engine.upper()) + '"'),
        envsh,
        flags=re.MULTILINE,
    )

    mode = get_choice(
        "mode",
        lambda: choice(
            "How will you run SubPEx?", choices=["MPI", "GPU", "MULTITHREAD"]
        ),
    )

    envsh = re.sub(
        r'^export MODE=".*?"',
        (('export MODE="' + mode.upper()) + '"'),
        envsh,
        flags=re.MULTILINE,
    )

    clear()
    return envsh


def define_adaptive_bins(adaptivepy, pcoord):
    # If pcoord is "jd", max is 1.0.
    log(f"STEP {get_step_num()}: Define adaptive bins", "-")
    log(
        "SubPEx uses WESTPA to enhance the sampling of your pocket. A progress coordinate assesses the extent of sampling. This progress coordinate is divided into bins, and WESTPA works to make sure the simulations collectively sample all bins equally."
    )

    log(f"You previously selected the progress coordinate: {pcoord}")
    if pcoord == "jd":
        maxcap = 1.0
        save_choice("maxcap", maxcap)
    else:
        log(
            "What is the maximum value of the progress coordinate beyond which SubPEx should no longer enhance sampling (in Angstroms)? We recommend 5."
        )

        maxcap = get_choice("maxcap", lambda: get_number("Maximum value"))
    log("\nHow many bins would you like to use? We recommend 15.")
    binsperdim = get_choice("binsperdim", lambda: int(get_number("Bin count")))
    log(
        "\nHow many walkers (mini simulations) would you like run per bin? Using more improves sampling at the cost of computer resources. We recommend 3."
    )

    bintargetcount = get_choice(
        "bintargetcount", lambda: int(get_number("Number of walkers per bin"))
    )

    adaptivepy = re.sub(
        r"^maxcap\s*=\s*\[.*\]",
        "maxcap = [" + str(maxcap) + "]",
        adaptivepy,
        flags=re.MULTILINE,
    )

    adaptivepy = re.sub(
        r"^binsperdim\s*=\s*\[.*\]",
        "binsperdim = [" + str(binsperdim) + "]",
        adaptivepy,
        flags=re.MULTILINE,
    )

    adaptivepy = re.sub(
        r"^bintargetcount\s*=\s*\d*",
        "bintargetcount = " + str(bintargetcount),
        adaptivepy,
        flags=re.MULTILINE,
    )

    clear()
    return adaptivepy


def update_run_time_and_job_name(westcfg):
    log(f"STEP {get_step_num()}: Update job run time and name", "-")
    log(
        "How long should your SubPEx job run? (Format your response like HH:MM:SS, e.g., 72:00:00 for 72 hours)"
    )
    run_time = get_choice("run_time", lambda: get_time("Run time"))

    log('\nWhat is the name of your SubPEx job? (e.g., "my_job")')
    job_name = get_choice("job_name", lambda: input("Job name: ").replace('"', ""))

    westcfg = re.sub(
        r"\bmax_run_wallclock:\s*\d{1,2}:\d{1,2}:\d{1,2}\b",
        "max_run_wallclock:    " + run_time,
        westcfg,
        flags=re.MULTILINE,
    )

    for submit_file in glob.glob("run.slurm*.sh"):
        with openfile(submit_file) as f:
            runsh = f.read()
        runsh = re.sub(
            r"^#SBATCH --time=.*$",
            "#SBATCH --time=" + run_time,
            runsh,
            flags=re.MULTILINE,
        )
        runsh = re.sub(
            r"^#SBATCH --job-name=.*$",
            "#SBATCH --job-name=" + job_name,
            runsh,
            flags=re.MULTILINE,
        )
        with open(submit_file, "w") as f:
            f.write(runsh)

    clear()
    return westcfg


def run_init():
    log(f"STEP {get_step_num()}: Initialize the SubPEx run", "-")
    if os.path.exists("west.h5"):
        log(
            "Would you like to clear the previous SubPEx run? WARNING: This will delete previously generated files, erasing the previous results so you can start a fresh SubPEx run."
        )

        if choice("Delete previous run and start over") == "y":
            print("")
            run_cmd(["rm", "-f", "./job_logs/*"])
            run_cmd(["rm", "-f", "west.h5"])
            run_cmd(["./init.sh"])
            if not os.path.exists("west.h5"):
                log("ERROR: Could not initialize the SubPEx run.")
                sys.exit(1)
            else:
                clear()
    else:
        import pdb

        pdb.set_trace()
        run_cmd(["rm", "-f", "./job_logs/*"])
        run_cmd(["./init.sh"])
        if not os.path.exists("west.h5"):
            log("ERROR: Could not initialize the SubPEx run.")
            sys.exit(1)
        else:
            clear()
    return westcfg


def finished():
    log("It appears the wizard was successful.")
    log("Optional steps NOT performed", "-")
    log(
        "The wizard has NOT performed any of the following optional steps. (The defaults should be fine in most cases.)"
    )
    log(
        '1. No changes to the list of auxiliary data to calculate ("west.cfg"). Unless you have modified "west.cfg", all such data will be calculated.'
    )
    log(
        '2. No changes to the "subpex -> resolution" field ("west.cfg"). If you plan to use the "jd" progress coordinate, you may wish to modify this field.'
    )
    log(
        '3. No changes to some of the minor parameters that control adaptive binning (e.g., "pcoordlength"; see "./adaptive_binning/adaptive.py").'
    )
    log("Critical steps NOT performed", "-")
    log(
        "The wizard has also NOT performed any of the following CRITICAL steps, which you must perform separately."
    )
    log(
        '1. Some needed changes still required to the "env.sh" file. Edit "env.sh" per your specific computing environment.'
    )
    log(
        '2. Some needed changes still required to the "run.slurm.*.sh" files. If using SLURM, edit these files per your environment.'
    )


# Load files
with openfile("./west.cfg") as f:
    westcfg = f.read()

with openfile("./adaptive_binning/adaptive.py") as f:
    adaptivepy = f.read()

with openfile("./westpa_scripts/get_pcoord.sh") as f:
    get_pcoord = f.read()

with openfile("./westpa_scripts/runseg.sh") as f:
    runseg = f.read()

with openfile("./env.sh") as f:
    envsh = f.read()

clear()
log(
    "This wizard is designed to help users setup and run SubPEx. See the README.md for more details."
)
log("You may exit this wizard at any time by pressing Ctrl+C.")

if check_if_restart_sim():
    # Exit
    sys.exit(0)

check_existing_params()
# confirm_environment()
# confirm_dependencies()
engine, get_pcoord, runseg = amber_or_namd(get_pcoord, runseg)
get_prelim_sim_files(engine)
get_sim_last_frame()
westcfg = update_westcfg_path_vars(westcfg)
westcfg, pcoord = pick_pcoord(westcfg)
westcfg = define_pocket(westcfg)
envsh = update_envsh(envsh, engine)
adaptivepy = define_adaptive_bins(adaptivepy, pcoord)
westcfg = update_run_time_and_job_name(westcfg)

# Save west.cfg
with open("./west.cfg", "w") as f:
    f.write(westcfg)

with open("./adaptive_binning/adaptive.py", "w") as f:
    f.write(adaptivepy)

with open("./westpa_scripts/get_pcoord.sh", "w") as f:
    f.write(get_pcoord)

with open("./westpa_scripts/runseg.sh", "w") as f:
    f.write(runseg)

with open("./env.sh", "w") as f:
    f.write(envsh)

run_init()
finished()
