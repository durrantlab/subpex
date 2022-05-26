# SubPEx 

## What is it?

Subpocket explorer (SubPEx) is a tool that uses weighted ensemble (WE), as implemented in 
WESTPA, to maximize pocket conformational search with the goal of obtaining an 
ensemble of protein conformations for their use in ensemble docking.

As any WE implementation it uses a progress coordinate to focus computational power on sampling phase space. The progress coordinates that we can use are:

- Jaccard Distance of pocket volumes (jd)
- backbone RMSD
- pocket heavy atoms RMSD
- composite RMSD (linear combination of backbone and pocket HA RMSD)

And we recommend using the composite RMSD.

## How to use it?

We have provided a yaml file to create your conda environment with all the dependencies that are needed. If you decide to make your own environment or you have westpa already installed SubPEx needs the following packages for the calculation of the progress coordinate: 

- MDAnalysis 
- NumPy
- SciPy
- scikit-learn
- yaml

At the moment manual input is needed to set-up the simualtion on the protein of interest. Hopefully this will be changed soon.

## How to run SubpEx

Before we start with the steps, we highly recommend visually inspecting the selected pocket. A way of visualizing the center of the pocket is to create a PDB file with one CA atom at the coordinates of your pocket and load it into your prefered visualization software (ChimeraX, PyMol, VMD, etc.) with the initial state of the protein of interest.

1. Clone or copy the repository (https://git.durrantlab.pitt.edu/erh91/SubPEx-Erich.git).
2. Soft link or copy the equilibrated trajectories and necessary restart files to the reference directory.
    - If using NAMD the dcd file is fine.
    - If using Amber the filetype that works with the SubPEx algorithm is the .nc.
3. Extract the last frame of the equilibrated trajectory.
4. Find the coordinates for center of the pocket and the radius you want to use.
    - you will be able to change this later once you visually inspect it.
5. Open the west.cfg file and modify it. 
    - add the center, radius and resolution.
    - select which progress coordinate and auxiliary data you want to calculate.
    - make sure that the WESTPA progress coordinate and auxdata matches the SubPEx ones.
    - make sure the paths to selection_file, topology, west_home, reference and reference_fop exist and are valid.
    - the reference is the pdb file that will used in EVERY SINGLE progress coordinate calculation.
    - selection_file and reference_fop will be calculated later but still need the names that will be used.
6. Run the westpa_scripts/get_reference_fop.py script. It uses the west.cfg as configuration file. Here is where the selection_file and reference_fop files are generated. 
7. Visually inspect the pocket field of points and the selection string (it uses MDAnalysis syntax) you can re-run the westpa_scripts/get_reference_fop.py script to recalculate them to fine tune your pocket.
8. Change the adaptive_binning/adaptive.py file to reflect the bins minimum and maximum values the progress coordinate must take and the number of bins and dimensions. 
9. Revise env.sh, you need to export the MD engine, this file can get complicated if using a supercomputing center. There are some extra files in the env-crc directory that will help with setup at a supercomputing center. We used the CRC (University of Pittsburgh's supercomputing center) and bridges2 of XSEDE and these files needed minimal changes.
10. Change westpa_scripts/get_pcoord.sh.
    - Make sure to link the files needed to start the simulations.
    - make sure that you are copying all the auxdata files to their respective WESTPA variable.
11. Activate the westpa conda environment and run the init.sh file.
12. If errors occur, check the job_logs directory. If there is no job_logs directory that alone causes it to fail.
13. Modify the runseg.sh
    - link files necessary files to restart simulation.
    - link all parent auxdata files as parent_AUXDATA.txt
    - make sure the random seed number gets added to the configuration file (comment out the one you are not using).
    - make sure the way to call the MD engine is correct.
    - make sure to call the right trajectory file for the jdistance.py script.
    - make sure that you are copying all the auxdata files to the westpa variables.
14. Modify the MD configuration file in reference (check names and other parameters)
    - Make sure the number of frames corresponds to pcoordlength minus one (pcoordlength is in the adaptive_binning/adaptive.py file)
15. run the run.sh file (unless you are running it in a supercomputing center).

__Notes__:
There are a LOT of parts that need to be perfect for it to run, since WESTPA sims are not that easy to set up. To debug check the west.log file and also the output in the job_logs directory.

## Important scripts and files and what they do 

- __env.sh__ sets up some environmental variables.
- __init.sh__ initializes the run. Creates the basis state, the initial states and 
calculates progress coordinates for those, using the bstate.py script.
- __gen_istate.sh__ makes istates directory.
- __get_pcoord.sh__ copies and links necessary files. Then calls on bstate.py script.
This script gets called by _init.sh_.
- __run.sh__ starts the run.
- __runseg.sh__ WESTPA runs this script for each trajectory segment. The script has 
three jobs:
    1) Link the necessary files for md simulations
    2) Run the MD simulation
    3) Calculate the progress coordinate using the jdistance.py script.
- __west.cfg__ the file that contains the configuration of the run.
- __west.h5__ contains all the results of the WE run.
