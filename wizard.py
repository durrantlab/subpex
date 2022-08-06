# The programs helps set up the SubPEx run. Runs only on linux.

import textwrap
import os
import json
import sys
import subprocess
import time

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
# YAML_LOADER = yaml.SafeLoader
YAML_LOADER = yaml.Loader

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
    log("Using previous choice: " + str(config[key]))
    return config[key]

def run_cmd(cmd):
    print("RUN: " + cmd + "\n")
    os.system(cmd)

def link_to_reference_mol(flnm):
    ext = flnm.split(".")[-1]
    run_cmd("ln -s " + flnm + " " + CUR_PATH + "/reference/mol." + ext)

def openfile(filename):
    if not os.path.exists(filename + ".orig"):
        os.system("cp " + filename + " " + filename + ".orig")

    return open(filename + ".orig", "r")

clear()

log("This wizard is designed to help users setup and run SubPEx. See the README.md for more details.")

log("You may exit this wizard at any time by pressing Ctrl+C.")

def confirm_environment():
    log("STEP 1: Environment", "-")
    log("You should run this wizard using an appropriate anaconda environment. To create the environment, install anaconda and create/activate the environment like this:")

    log("conda env create -f environment.yaml\n conda activate westpa")

    if choice("Have you created and activated the environment?") == "n":
        log("\nPlease create and activate the environment. Then try again.")
        exit(1)

    clear()

def confirm_dependencies():
    log("Step 2: Check dependencies", "-")

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

def amber_or_namd():
    log("STEP 3: MD engine", "-")
    engine = get_choice(
        "engine", 
        lambda : choice(
            "Will you use the Amber or NAMD MD engine?", 
            choices=["Amber", "NAMD"]
        )
    )
    clear()
    return engine

def get_prelim_sim_files(engine):
    log("STEP 4: Preliminary MD files", "-")
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
    run_cmd("rm -f " + CUR_PATH + "/reference/mol.*")

    # Link files
    link_to_reference_mol(coor_file)
    for restart_file in restart_files:
        link_to_reference_mol(restart_file)
    link_to_reference_mol(prmtop_file)
    clear()
    
def get_sim_last_frame():
    log("STEP 6: Preliminary MD last frame", "-")
    log("SubPEx needs the last frame of the preliminary trajectory as a `pdb` file. You can use programs such as VMD to save the last frame.")

    last_frame_file = get_choice(
        "last_frame_file",
        lambda : file_path("Last-frame .pdb file", ".pdb")
    )

    link_to_reference_mol(last_frame_file)
    clear()

def update_westcfg_path_vars():
    # Open west.cfg yaml file
    f = openfile("west.cfg")
    westcfg = yaml.load(f, Loader=YAML_LOADER)
    westcfg["subpex"]["selection_file"] = CUR_PATH + "/reference/selection_string.txt"
    westcfg["subpex"]["topology"] = CUR_PATH + "/reference/mol.prmtop"
    westcfg["subpex"]["west_home"] = CUR_PATH
    westcfg["subpex"]["reference"] = CUR_PATH + "/reference/last_frame.pdb"
    westcfg["subpex"]["reference_fop"] = CUR_PATH + "/reference/fop_ref.xyz"

    return westcfg

def pick_pcoord(westcfg):
    log("STEP 7: Progress coordinate", "-")
    log("Which progress coordinate would you like to use?")
    log("composite: combination of pocket and back-bone RMSD (recommended)\nprmsd: RMSD of pocket-lining atoms\nbb: RMSD of all back-bone atoms\njd: Jaccard distance between pocket shapes")

    pcoord = get_choice(
        "pcoord", 
        lambda : choice(
            "Which progress coordinate?", 
            choices=["composite", "prmsd", "bb", "jd"]
        )
    )
    westcfg["subpex"]["pcoord"] = pcoord
    clear()

def pick_auxdata(westcfg):
    log("STEP 8: Select auxiliary data", "-")
    print("STILL NEED TO IMPLEMENT!!!")  # TODO:
    clear()

def define_pocket(westcfg):
    log("STEP 9: Define the binding pocket", "-")
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

    westcfg["subpex"]["center"] = [x_coor, y_coor, z_coor]
    westcfg["subpex"]["radius"] = pocket_radius

    # Save west.cfg
    with open("west.cfg", "w") as f:
        yaml.dump(westcfg, f)

    # Generate xyz file
    prts = [PYTHON_EXE, "./westpa_scripts/get_reference_fop.py", "west.cfg"]
    # prts = ["ls"]
    cmd = subprocess.Popen(
        prts, 
        stdin=subprocess.PIPE, stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE, text=True
    )
    while cmd.poll() is None:
        time.sleep(0.1)
        output, errors = cmd.communicate()
        print(output)
    cmd.wait()
    print(output, errors)

    # os.system(PYTHON_EXE + " ./westpa_scripts/get_reference_fop.py west.cfg")

    # Note that intentionally not allowing user to specify FOP resolution here.
    # Keep it simple.
    # clear()
    

# confirm_environment()
# confirm_dependencies()
engine = amber_or_namd()
get_prelim_sim_files(engine)
get_sim_last_frame()
westcfg = update_westcfg_path_vars()
pick_pcoord(westcfg)
pick_auxdata(westcfg)
define_pocket(westcfg)

# Save west.cfg
with open("west.cfg", "w") as f:
    yaml.dump(westcfg, f)