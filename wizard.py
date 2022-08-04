# The programs helps set up the SubPEx run. Runs only on linux.

import textwrap
import os
import json

try:
    import yaml
except:
    print("\nPyYAML is missing. Have you set up your westpa conda environment?\n")
    print("conda env create -f environment.yaml\nconda activate westpa\n")
    exit(1)

if os.name == 'nt':
    print("Wizard script supported only on Linux. Sorry!")
    exit(1)

CUR_PATH = os.path.abspath(os.path.dirname(__file__))

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
    log("Using previous choice: " + config[key])
    return config[key]

def run_cmd(cmd):
    print("RUN: " + cmd + "\n")
    os.system(cmd)

def link_to_reference_mol(flnm):
    ext = flnm.split(".")[-1]
    run_cmd("ln -s " + flnm + " " + CUR_PATH + "/reference/mol." + ext)

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
        env = yaml.load(f, Loader=yaml.SafeLoader)
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

confirm_environment()
confirm_dependencies()
engine = amber_or_namd()
get_prelim_sim_files(engine)
get_sim_last_frame()
