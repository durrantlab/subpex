# The programs helps set up the SubPEx run. Runs only on linux.

import textwrap
import os
import json
import sys
import subprocess
import time
import re

try:
    import yaml
except:
    print("\nPyYAML is missing. Have you set up your westpa conda environment?\n")
    print("conda env create -f environment.yaml\nconda activate westpa\n")
    exit(1)

if os.name == 'nt':
    print("\nWizard script supported only on Linux. Sorry!\n")
    exit(1)

# Make sure west.cfg in present directory
if not os.path.exists("west.cfg"):
    print("\nPlease run this script from the same directory containing the west.cfg file.\n")
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

def choice(prompt, choices=["y", "n"]):
    lowers = [c.lower() for c in choices]
    while True:
        answer = input(prompt + " (" + "/".join(choices) + "): ").lower()
        if answer not in lowers:
            if len(choices) == 2:
                print("Please answer " + " or ".join(choices))
            else:
                print("Please answer " + ", ".join(choices[:-1]) + ", or " + choices[-1])
            continue
        return answer

def clear():
    os.system('clear')

def file_path(prompt, ext):
    while True:
        answer = input(prompt + ": ").strip()
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
        answer = input(prompt + ": ").strip()

        # Check if it's a number
        try:
            answer = float(answer)
        except:
            log("\nPlease enter a number. Try again.")
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
    log("Using previous choice for " + key + ": " + str(config[key]))
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

def run_cmd(prts):
    log("Running command: " + " ".join(prts))
    cmd = subprocess.Popen(
        prts, 
        stdin=subprocess.PIPE, stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE, text=True
    )
    cmd.wait()
    output, errors = cmd.communicate()
    output = output.strip()
    errors = errors.strip()

    if output:
        log(output)
    if errors:
        log(errors)

def link_to_reference_mol(flnm):
    ext = flnm.split(".")[-1]
    run_cmd(["ln", "-s", flnm,  CUR_PATH + "/reference/mol." + ext])

def openfile(filename):
    if not os.path.exists(filename + ".orig"):
        os.system("cp " + filename + " " + filename + ".orig")

    return open(filename + ".orig", "r")

def confirm_environment():
    log("STEP " + get_step_num() + ": Environment", "-")
    log("You should run this wizard using an appropriate anaconda environment. To create the environment, install anaconda and create/activate the environment like this:")

    log("conda env create -f environment.yaml\n conda activate westpa")

    if choice("Have you created and activated the environment?") == "n":
        log("\nPlease create and activate the environment. Then try again.")
        exit(1)

    clear()

def confirm_dependencies():
    log("Step " + get_step_num() + ": Check dependencies", "-")

    # Load environment.yaml into a dictionary
    with open("environment.yaml", "r") as f:
        env = yaml.load(f, Loader=YAML_LOADER)
    env = [v for v in env["dependencies"] if v[:1] != "_"]
    env = [v.split("=")[0] for v in env]
    to_ignore=["biopython", "bzip2", "curl", "freetype", "griddataformats", "hdf4", "hdf5", "ipython", "jpeg", "krb5", "lcms2", "libblas", "libcblas", "libcurl", "libedit", "libev", "libffi", "libgfortran5", "libgomp", "liblapack", "libnetcdf", "libnghttp2", "libopenblas", "libpng", "libsodium", "libssh2", "libtiff", "mdanalysis", "ncurses", "netcdf4", "openssl", "pillow", "python", "python_abi", "pyyaml", "pyzmq", "sqlite", "tk", "westpa", "xz", "zeromq", "zstd"]

    for v in env:
        if v in to_ignore:
            continue
        if "-" in v:
            continue
        print("Testing module " + v)
        # import by name
        try:
            __import__(v)
        except:
            log("\nModule " + v + " is missing. Have you created and activated the westpa environment?")
            log("conda env create -f environment.yaml\n conda activate westpa")
            log("Please create/activate it and try again.")
            exit(1)
    clear()

def amber_or_namd(get_pcoord, runseg):
    log("STEP " + get_step_num() + ": MD engine", "-")
    engine = get_choice(
        "engine", 
        lambda : choice(
            "Will you use the Amber or NAMD MD engine?", 
            choices=["Amber", "NAMD"]
        )
    )
    clear()

    get_pcoord = re.sub(r'^export ENGINE=".*?"', 'export ENGINE="' + engine.upper() + '"', get_pcoord, flags=re.MULTILINE)
    runseg = re.sub(r'^export ENGINE=".*?"', 'export ENGINE="' + engine.upper() + '"', get_pcoord, flags=re.MULTILINE)

    return engine, get_pcoord, runseg

def get_prelim_sim_files(engine):
    log("STEP " + get_step_num() + ": Preliminary MD files", "-")
    log("SubPEx assumes you have already run preliminary simulations to equilibrate your system. Where are the coordinate and restart files associated with these preliminary simulations?")
    log("(Note that you can press Ctrl-Z to pause the wizard and search for the file. Type `fg` to resume the wizard when done.)")
    
    if engine == "amber":
        coor_file = get_choice(
            "coor_file", 
            lambda : file_path("Amber .nc file", ".nc")
        )
        restart_files = [get_choice(
            "restart_files", 
            lambda : file_path("Amber .rst file", ".rst")
        )]
    elif engine == "namd":
        coor_file = get_choice(
            "coor_file", 
            lambda : file_path("NAMD .dcd file", ".dcd")
        )
        restart_files = []
        for ext in [".xsc", ".coor", ".inpcrd"]:
            restart_files.append(
                get_choice(
                    "_ignore", 
                    lambda : file_path("NAMD " + ext + " file", ext)
                )
            )
        save_choice("restart_files", restart_files)

    # Both engines require a .prmtop file
    prmtop_file = get_choice(
        "prmtop_file",
        lambda : file_path("Parameter (.prmtop) file", ".prmtop")
    )

    # Remove any old mol files
    run_cmd(["rm", "-f", CUR_PATH + "/reference/mol.*"])

    # Link files
    link_to_reference_mol(coor_file)
    for restart_file in restart_files:
        link_to_reference_mol(restart_file)
    link_to_reference_mol(prmtop_file)
    clear()
    
def get_sim_last_frame():
    log("STEP " + get_step_num() + ": Preliminary MD last frame", "-")
    log("SubPEx needs the last frame of the preliminary trajectory as a `pdb` file. You can use programs such as VMD to save the last frame.")

    last_frame_file = get_choice(
        "last_frame_file",
        lambda : file_path("Last-frame .pdb file", ".pdb")
    )

    link_to_reference_mol(last_frame_file)
    clear()

def update_westcfg_path_vars(westcfg):
    # Open west.cfg yaml file
    westcfg = re.sub(r"\bselection_file: .*?$", "selection_file: " + CUR_PATH + "/reference/selection_string.txt", westcfg, flags=re.MULTILINE)
    westcfg = re.sub(r"\btopology: .*?$", "topology: " + CUR_PATH + "/reference/mol.prmtop", westcfg, flags=re.MULTILINE)
    westcfg = re.sub(r"\bwest_home: .*?$", "west_home: " + CUR_PATH, westcfg, flags=re.MULTILINE)
    westcfg = re.sub(r"\breference: .*?$", "reference: " + CUR_PATH + "/reference/mol.pdb", westcfg, flags=re.MULTILINE)
    westcfg = re.sub(r"\breference_fop: .*?$", "reference_fop: " + CUR_PATH + "/reference/fop_ref.xyz", westcfg, flags=re.MULTILINE)

    return westcfg

def pick_pcoord(westcfg):
    log("STEP " + get_step_num() + ": Progress coordinate", "-")
    log("Which progress coordinate would you like to use?")
    log("composite: combination of pocket and back-bone RMSD (recommended)\nprmsd: RMSD of pocket-lining atoms\nbb: RMSD of all back-bone atoms\njd: Jaccard distance between pocket shapes")

    pcoord = get_choice(
        "pcoord", 
        lambda : choice(
            "Which progress coordinate?", 
            choices=["composite", "prmsd", "bb", "jd"]
        )
    )
    westcfg = re.sub(r"(\bpcoord:\s*\n\s*- ).*?$", r"\1" + pcoord, westcfg, flags=re.MULTILINE)
    clear()
    return westcfg, pcoord

def pick_auxdata(westcfg):
    log("STEP " + get_step_num() + ": Select auxiliary data", "-")
    print("STILL NEED TO IMPLEMENT!!!")  # TODO:
    clear()
    return westcfg

def define_pocket(westcfg):
    log("STEP " + get_step_num() + ": Define the binding pocket", "-")
    log("You must define the location of the binding pocket you wish to sample. Visual inspection (e.g., using VMD) is often useful for this step. See README.md for more suggestions.")
    log("Be sure to specify the location of the pocket relative to the reference PDB file: " + CUR_PATH + "/reference/mol.pdb")
    
    pocket_radius = get_choice(
        "pocket_radius", 
        lambda : get_number("Pocket radius")
    )

    x_coor = get_choice(
        "x_coor", 
        lambda : get_number("X coordinate of pocket center")
    )

    y_coor = get_choice(
        "y_coor", 
        lambda : get_number("Y coordinate of pocket center")
    )

    z_coor = get_choice(
        "z_coor", 
        lambda : get_number("Z coordinate of pocket center")
    )

    westcfg = re.sub(r"\bcenter: \[.*?\]", "center: [" + str(x_coor) + ", " + str(y_coor) + ", " + str(z_coor) + "]", westcfg, flags=re.MULTILINE)
    westcfg = re.sub(r"\bradius: .*?$", "radius: " + str(pocket_radius), westcfg, flags=re.MULTILINE)

    # Save west.cfg
    with open("west.cfg", "w") as f:
        f.write(westcfg)

    # Generate xyz file
    run_cmd(["rm", "-f", CUR_PATH + "/reference/fop_ref.xyz"])
    run_cmd(["rm", "-f", CUR_PATH + "/reference/selection_string.txt"])
    run_cmd([PYTHON_EXE, "./westpa_scripts/get_reference_fop.py", "west.cfg"])

    # Check if file exists
    if not os.path.exists(CUR_PATH + "/reference/fop_ref.xyz") or not os.path.exists(CUR_PATH + "/reference/selection_string.txt"):
        log("ERROR: Could not generate reference-fop file (fop_ref.xyz) and/or selection-string file (selection_string.txt) in the directory " + CUR_PATH + "/reference/")
        log("Are you sure the center and radius specified encompass at least some of the pocket-lining residues?")
        sys.exit(1)

    # Note that intentionally not allowing user to specify FOP resolution here.
    # Keep it simple.
    clear()

    westcfg = check_pocket(westcfg)

    return westcfg

def check_pocket(westcfg):
    log("STEP " + get_step_num() + ": Check the binding pocket", "-")
    
    log("The wizard has generated a pocket-field-of-points (fop) file (" + CUR_PATH + "/reference/fop_ref.xyz) and a selection-string file (" + CUR_PATH + "/reference/selection_string.txt).")
    log("Please visually inspect these files to:")
    log("1. Ensure that the points in fop_ref.xyz entirely fill the pocket of interest.")
    log("2. Ensure that the residues described in selection_string.txt truly line the pocket of interest.")
    log("(Note that the popular molecular visualization program VMD can load xyz files and select residues.")
    
    if choice("Can you confirm the pocket is properly defined?") == "n":
        clear()
        forget_choice("pocket_radius")
        forget_choice("x_coor")
        forget_choice("y_coor")
        forget_choice("z_coor")
        westcfg = define_pocket(westcfg)

    clear()
    return westcfg

def define_adaptive_bins(adaptivepy, pcoord):
    # If pcoord is "jd", max is 1.0.
    log("STEP " + get_step_num() + ": Define adaptive bins", "-")
    
    log("SubPEx uses WESTPA to enhance the sampling of your pocket. A progress coordinate assesses the extent of sampling. This progress coordinate is divided into bins, and WESTPA works to make sure the simulations collectively sample all bins equally.")

    log("You previously selected the progress coordinate: " + pcoord)

    if pcoord == "jd":
        maxcap = 1.0
        save_choice("maxcap", maxcap)
    else:
        log("What is the maximum value of the progress coordinate beyond which SubPEx should no longer enhance sampling (in Angstroms)? We recommend 5.")
        maxcap = get_choice(
            "maxcap", 
            lambda : get_number("Maximum value")
        )

    log("\nHow many bins would you like to use? We recommend 15.")
    binsperdim = get_choice(
        "binsperdim", 
        lambda : int(get_number("Bin count"))
    )

    log("\nHow many walkers (mini simulations) would you like run per bin? Using more improves sampling at the cost of computer resources. We recommend 3.")
    bintargetcount = get_choice(
        "bintargetcount", 
        lambda : int(get_number("Number of walkers per bin"))
    )

    adaptivepy = re.sub(r"^maxcap\s*=\s*\[.*?\]", "maxcap = [" + str(maxcap) + "]", adaptivepy, flags=re.MULTILINE)
    adaptivepy = re.sub(r"^binsperdim\s*=\s*\[.*?\]", "binsperdim = [" + str(binsperdim) + "]", adaptivepy, flags=re.MULTILINE)
    adaptivepy = re.sub(r"^bintargetcount\s*=\s*[0-9]*", "bintargetcount = " + str(bintargetcount), adaptivepy, flags=re.MULTILINE)

    clear()
    return adaptivepy

def run_init():
    log("STEP " + get_step_num() + ": Initialize the SubPEx run", "-")

    if os.path.exists("west.h5"):
        log("Would you like to clear the previous SubPEx run? WARNING: This will delete previously generated files, erasing the previous results so you can start a fresh SubPEx run.")

        if choice("Delete previous run and start over") == "y":
            print("")
            run_cmd(["rm", "-f", "west.h5"])
            run_cmd(["./init.sh"])
            if not os.path.exists("west.h5"):
                log("ERROR: Could not initialize the SubPEx run.")
                sys.exit(1)
            else:
                clear()
    else:
        run_cmd(["./init.sh"])
        if not os.path.exists("west.h5"):
            log("ERROR: Could not initialize the SubPEx run.")
            sys.exit(1)
        else:
            clear()

    return westcfg

# Load files
with openfile("./west.cfg") as f:
    westcfg = f.read()

with openfile("./adaptive_binning/adaptive.py") as f:
    adaptivepy = f.read()

with openfile("./westpa_scripts/get_pcoord.sh") as f:
    get_pcoord = f.read()

with openfile("./westpa_scripts/runseg.sh") as f:
    runseg = f.read()

clear()
log("This wizard is designed to help users setup and run SubPEx. See the README.md for more details.")
log("You may exit this wizard at any time by pressing Ctrl+C.")

# confirm_environment()
# confirm_dependencies()
engine, get_pcoord, runseg = amber_or_namd(get_pcoord, runseg)
get_prelim_sim_files(engine)
get_sim_last_frame()
westcfg = update_westcfg_path_vars(westcfg)
westcfg, pcoord = pick_pcoord(westcfg)
westcfg = pick_auxdata(westcfg)
westcfg = define_pocket(westcfg)
adaptivepy = define_adaptive_bins(adaptivepy, pcoord)
run_init()

# appropriate WORKMANAGER (env.sh)
# Summary of next steps, and customizations not made.

# Save west.cfg
with open("./west.cfg", "w") as f:
    f.write(westcfg)

with open("./adaptive_binning/adaptive.py", "w") as f:
    f.write(adaptivepy)

with open("./westpa_scripts/get_pcoord.sh", "w") as f:
    f.write(get_pcoord)

with open("./westpa_scripts/runseg.sh", "w") as f:
    f.write(runseg)