# SubPEx 

## What is it?

Subpocket explorer (SubPEx) is a tool that uses weighted ensemble (WE), as
implemented in WESTPA, to maximize pocket conformational search with the goal of
obtaining an ensemble of protein conformations for use in ensemble docking.

As with any WE implementation, SubPEx uses a progress coordinate to focus
computational power on sampling phase space. The available progress coordinates
are:

- Jaccard distance of pocket volumes (jd)
- backbone RMSD
- pocket heavy atoms RMSD
- composite RMSD (a linear combination of backbone and pocket heavy-atom RMSD)

We recommend using the composite RMSD progress coordinate.

## How to use run SubPEx?

The first step to start running SubPEx simulations is to download this repository, this can be achieved by cloning or copying the repository. 

```bash
git clone https://git.durrantlab.pitt.edu/erh91/SubPEx-Erich.git
```

In the reposirtory we have provided a yaml file to create a `conda` environment with all the necessary
dependencies. Simply run this command:

```bash
conda env create -f environment.yaml
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

Manual input is needed to set up the simulation of the protein of interest.
Hopefully, this will be changed soon and an autobuilder will be available.

1. Soft link or copy the equilibrated trajectories and necessary restart files
   to the reference directory.
    - If using NAMD, the .dcd file is fine.

    - If using Amber, the filetype that works with the SubPEx algorithm is the
      .nc.
```bash
ln -s \file\path\to\trajectory\files \WEST\ROOT\reference
```
2. Extract the last frame of the equilibrated trajectory with your prefered package.
3. Find the coordinates for the center of the pocket and the radius you want to
   use.
    - you will be able to change this later once you visually inspect it.
    - visual inspection can be done by creating a PDB file with the with one
      CA atom at the coordinates of your pocket and load it into your preferred
      visualization software (ChimeraX, PyMol, VMD, etc.) with the extracted last frame 
      of the previous step.
4. Open the west.cfg file and modify it. 
    - add the center, radius, and resolution.
    - select which progress coordinate and auxiliary data you want to calculate.
    - make sure that the WESTPA progress coordinate and auxdata match the SubPEx
      ones.
    - make sure the paths to selection_file, topology, west_home, reference, and
      reference_fop exist and are valid.
    - the reference is the PDB file that will be used in EVERY SINGLE progress
      coordinate calculation (the one obtained in step 2).
    - selection_file and reference_fop will be calculated later but still need
      the names used.
5. Run the westpa_scripts/get_reference_fop.py script. It uses the west.cfg as
   the configuration file. Here is where the selection_file and reference_fop
   files are generated. 
6. Visually inspect the pocket field of points and/or the selection string (it uses
   MDAnalysis syntax). You can re-run the westpa_scripts/get_reference_fop.py
   script to recalculate them to fine-tune your pocket.
7. Change the adaptive_binning/adaptive.py file to reflect the bins' minimum and
   maximum values the progress coordinate must take and the number of bins and
   dimensions. __NOTE__: We do not need target states, so we can ignore those parameters.
8. Revise env.sh. You need to export the MD engine, and this file can get
   complicated if using a supercomputing center. Some extra files in the env-crc
   directory will help with the setup at a supercomputing center. We used the
   CRC (University of Pittsburgh's supercomputing center) and bridges2 of XSEDE,
   and these files needed minimal changes.
9. Change westpa_scripts/get_pcoord.sh.
    - Make sure to link the files needed to start the simulations. (check lines 23-25)
      The file has example commands for a NAMD run.
    - make sure you are copying all the auxdata files to their respective WESTPA
      variable. (check lines 32-36)
10. Activate the WESTPA conda environment and run the init.sh file.
```bash
conda activate westpa
```
11. If errors occur, check the job_logs directory. If there is no job_logs
    directory, that alone causes it to fail.
12. Modify the runseg.sh
    - link necessary files to restart the simulation. (check lines 79-88)
    - link all parent auxdata files as parent_AUXDATA.txt (check lines 59-64)
    - make sure the random seed number gets added to the configuration file
      (comment out the one you are not using). (check lines 66-72)
    - make sure the way to call the MD engine is correct. (check lines 92-100)
    - make sure to call the correct trajectory file for the pcoord.py script. (check lines 110-112)
    - make sure you copy all the auxdata files to the WESTPA variables. (check lines 114-119)
13. Modify the MD configuration file in reference (check names and other
    parameters)
    - Make sure the number of frames corresponds to pcoordlength minus one
      (pcoordlength is in the adaptive_binning/adaptive.py file)
14. Run the run.sh file (unless you are running it in a supercomputing center).
    We have provided a template file (subpex.sh) that works for our supercomputing center, 
    it may not work for your center. Please check with your IT person to troubleshoot any
    problem.

__Notes__: A lot of parts need to be perfect for it to run since WESTPA sims are
not that easy to set up. To debug, check the west.log file and the output in the
job_logs directory.

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
