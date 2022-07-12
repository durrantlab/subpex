# SubPEx

## What is it?

Subpocket explorer (SubPEx) is a tool that uses weighted ensemble (WE), as
implemented in WESTPA, to maximize pocket conformational search. SubPEx's goal
is to produce a diverse ensemble of protein conformations for use in ensemble
docking.

As with any WE implementation, SubPEx uses a progress coordinate to focus
computational power on sampling phase space. The available progress coordinates
are:

- composite RMSD (a linear combination of backbone and pocket heavy-atom RMSD)
- pocket heavy atoms RMSD
- backbone RMSD
- Jaccard distance of pocket volumes (jd)

We recommend using the composite RMSD progress coordinate.

## SubPEx use

### Downloading and installing SubPEx

The first stepis to download, clone, or copy the repository.

```bash
git clone https://git.durrantlab.pitt.edu/erh91/SubPEx-Erich.git
```

The repository includes a `yaml` file to create a `conda` environment with all
the necessary dependencies. Simply run this command:

```bash
conda env create -f environment.yaml
```

To activate the new `conda` environment, run:

```bash
conda activate westpa
```

Some users may wish to create their own environments or to use an existing
WESTPA environment. If so, install the following packages so SubPEx can
calculate the progress coordinate:

- Westpa
- MDAnalysis
- NumPy
- SciPy
- scikit-learn
- yaml

### Configuring a SubPEx run

Currently, users must manually set up their SubPEx simulations by editing key
SubPEx/WESTPA files. We hope to soon develop an autobuilder ("wizard") to
simplify this process.

___Link your preliminary, equilibrated simulation___

1. SubPEx assumes you have already run preliminary simulations to equilibrate
   your system. Soft link or copy your preliminary, equilibrated trajectories
   and necessary restart files (***JDD: Possible to be more descriptive?***) to
   the `./reference/` directory. (Note that this directory already contains the
   `md.conf` and `prod_npt.in` template files, which SubPEX uses to interface
   with the NAMD and AMBER MD engines, respectively.)
    - If using NAMD, copying the `.dcd` file is fine.
    - If using Amber, the filetype that works with the SubPEx algorithm is the
      `.nc`.

      ```bash
      ln -s \file\path\to\trajectory\files \WEST\ROOT\reference\
      ```

2. Extract the last frame of the preliminary, equilibrated trajectory as a `pdb`
   file with your preferred molecular analysis program (e.g., VMD).

___Edit the `west.cfg` file___

1. Edit the following parameters in the `west.cfg` file:
    - the path variables:
       - `reference`: the PDB file that will be used in EVERY SINGLE
         progress-coordinate calculation (the last frame of the preliminary,
         equilibrated simulation mentioned above).
       - `selection_file`: path to a text file containing the pocket selection
         string (MDAnalysis selection notation). This file will be automatically
         generated in a subsequent step.
       - `reference_fop`: path to an `xyz` file containing the field of points
         needed to calculate the `jd` progress coordinate. This file is also
         useful for visualizing the selected pocket. It will be automatically
         generated in a subsequent step.
       - `west_home`: home directory of the SubPEx run. You'll most likely want
         to use the same directory that contains the `west.cfg` file itself.
       - `topology`: topology file needed for the MD simulations (likely the
         same topology file used in the preliminary, equilibrated simulations).
    - the progress coordinate (`pcoord`) to use.
       - `composite`: composite RMSD (recommended)
       - `prmsd`: pocket heavy atoms RMSD
       - `bb`: backbone RMSD
       - `jd`: Jaccard distance
    - the auxiliary data (`auxdata`) to calculate and save.
       - `composite`: composite RMSD
       - `prmsd`: pocket heavy atoms RMSD*
       - `pvol`: pocket volume
       - `bb`: backbone RMSD
       - `rog_pocket`: radius of gyration of the pocket
       - `jd`: Jaccard distance
    - make sure that the WESTPA progress coordinate and auxdata match the SubPEx
      ones (these sections are both found in the `west.cfg` file).
      - ***JDD: Where exactly? data -> datasets and also executable -> datasets?***

___Define the pocket to sample___

1. You must define the location of the binding pocket you wish to sample. Find
   the coordinates of the pocket center and radius using the extracted last
   frame.
    - Visual inspection is often useful at this step. You might create a PDB
      file with a CA dummy atom. Load that together with the extracted last
      frame of the previous step into your preferred visualization software
      (ChimeraX, PyMol, VMD, etc.).
    - Manually move the dummy atom to the pocket center and measure its
      location. Similarly, use the dummy atom to determine the radius from that
      center required to encompass the pocket of interest.
4. Return to the `west.cfg` file and edit the following parameters:
   - `center`: the pocket center
   - `radius`: the pocket radius
   - `resolution`: the distance between adjacent pocket-filling grid points
     (especially important if using the `jd` progress coordinate)
5. Run `python westpa_scripts/get_reference_fop.py west.cfg`. This script will
   generate the files specified by the `selection_file` and `reference_fop`
   parameters in the `west.cfg` file.
6. Visually inspect the pocket field of points (fop) and/or the selection string
   (MDAnalysis selection syntax).
   - Ensure that the fop (`reference_fop`) entirely fills the pocket of
     interest.
   - Ensure that the residues (`selection_file`) truly line the pocket of
     interest.
   - Note that the popular molecular visualization program VMD can load `xyz`
     files and select residues.
7. After visual inspection, adjust the `west.cfg` file (`center`, `radius`, and
     `resolution` parameters) and re-run the
     `westpa_scripts/get_reference_fop.py` script. Continue to recalculate the
     pocket as needed to fine-tune your pocket.

___Setup the progress coordinate calculations___

1. Update the variables in the `adaptive_binning/adaptive.py` file to indicate
   the number of walkers per bin, the bins' minimum and maximum values, etc.
    - This file controls the adaptive binning scheme that SubPEx uses.
    - A detailed description of each variable is given in the file itself.
2. Change `westpa_scripts/get_pcoord.sh`.
    - This script runs when calculating initial progress coordinates for new
      initial states (istates).
    - Make sure to link the files needed to start the SubPEx simulations (the
      ones you put in the `./reference/` directory above).
      - Check lines 24-26. Replace {REFERENCE}, {RESTART_FILE}, and
        {TOPOLOGY_FILE} with the appropriate file names.
      - The script includes examples for NAMD and AMBER runs.
    <!-- - make sure you are copying all the auxdata files to their respective WESTPA -->
      <!-- variable. (check lines 38-42) TODO: Hoping this isn't necessary. See relevant file. -->
12. Modify the `runseg.sh`
    - This file runs each WESTPA/SubPEx walker (segment). It creates the needed
      directory, runs the walker simulation, and calculates the progress
      coordinate.
    <!-- - link all parent auxdata files as parent_AUXDATA.txt (check lines 59-64) TODO: Hoping this isn't needed anymore. -->
    - Make sure the random seed number gets added to the MD-engine configuration
      file (check lines 66-72). Be sure to comment out the NAMD code and
      uncomment the AMBER code if using AMBER as the MD engine.
    - Link the necessary files to restart the walker simulation (check lines
      79-88). Be sure to comment out/in the appropriate block for NAMD vs.
      AMBER.
    - Be sure to call the correct MD engine (check lines 92-100).
    - Be sure to call the correct trajectory file for the pcoord.py script.
      (check lines 110-114)
    <!-- - Make sure you copy all the auxdata files to the WESTPA variables. (check
      lines 114-119) TODO: Hopefully not needed now. -->

8. Revise `env.sh`. You need to export the MD engine, and this file can get
   complicated if using a supercomputing center. Some extra files in the
   `./env-crc/` directory will help with the setup at a supercomputing center.
   We used the CRC (University of Pittsburgh's supercomputing center) and
   bridges2 of XSEDE, and these files needed minimal changes. ***JDD: Not clear
   what to change in this file. LOOK AT `./env-crc` directory and give thought
   to how to make more usable.***
10. Activate the WESTPA conda environment and run the init.sh file.

```bash
conda activate westpa
```

11. If errors occur, check the `./job_logs` directory. If there is no
    `./job_logs` directory, that alone causes it to fail.
13. Modify the MD configuration file in reference (check names and other
    parameters)
    - Make sure the number of frames corresponds to pcoordlength minus one
      (pcoordlength is in the `adaptive_binning/adaptive.py` file)
14. Run the run.sh file (unless you are running it in a supercomputing center).
    We have provided a template file (subpex.sh) that works for our
    supercomputing center, it may not work for your center. Please check with
    your IT person to troubleshoot any problem.

__Notes__: A lot of parts need to be perfect for it to run since WESTPA sims are
not that easy to set up. To debug, check the `west.log` file and the output in
the `job_logs` directory.

## Important scripts and files and what they do

- __env.sh__ sets up some environmental variables.
- __init.sh__ initializes the run. Creates the basis and initial states and
  calculates progress coordinates for those, using the bstate.py script.
- __gen_istate.sh__ makes istates directory.
- __get_pcoord.sh__ copies and links necessary files. Then calls on bstate.py
  script. This script gets called by _init.sh_.
- __run.sh__ starts the run.
- __runseg.sh__ WESTPA runs this script for each trajectory segment. The script
  has three jobs:
    1) Link the necessary files for md simulations
    2) Run the MD simulation
    3) Calculate the progress coordinate using the jdistance.py script.
- __west.cfg__ the file containing the run's configuration.
- __west.h5__ contains all the results of the WE run.
- __get_reference_fop.py__ calcualtes the initial field of points for JD progress 
  coordinate and creates selection string for MDAnalysis.
- __pcoord_istate.py__ calcualtes progres coordinate for the initial states.
- __pcoord.py__ calculates the progress coordiante for the production run.
