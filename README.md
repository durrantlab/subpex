# SubPEx code

## What is it?

Subpocket explorer (SubPEx) is a tool that uses weighted ensemble, as implemented in 
WESTPA, to maximizes pocket conformational search with the goal of obtaining an 
ensemble of protein conformations for the use in ensemble docking.

## How to use it?

Make sure you have WESTPA installed and that it calls the right python environment, or install the westpa conda environment. For the calculation of the progress coordinate SubPEx needs the following packages: 

- MDAnalysis 
- NumPy
- SciPy
- scikit-learn
- yaml or json

At the moment manual input is needed for the user to be able to adapt it to the protein
of interest. Hopefully this will be changed soon.

Before running anything you need to create the reference field of points and the selection string. To do so, you need to run the __westpa\_scripts/get\_reference\_fop.py__.

Then, run init.sh which will calculate the basis state progress coordinates and then create the initial states and 
calculate their progress coordinate. After that you will need to run the run.sh script.

For more information on these scripts go to the comments in the files.

## how to run SubpEx

1. Clone the repository (https://git.durrantlab.pitt.edu/erh91/SubPEx-Erich.git).
2. Soft link or copy the equilibrated trajectories and necessary restart files to the reference directory.
    - If using NAMD the dcd file is fine.
    - If using Amber the filetype that works with the SubPEx algorithm is the .nc.
3. Extract the last frame of the equilibrated trajectory.
4. Find the coordinates for center of the pocket and the radius you want to use.
    - you will be able to change this later once you visually inspect it.
5. Open the west.cfg file and modify it.
    - add the center, radius and resolution.
    - Make sure that the WESTPA progress coordinate and auxdata matches the SubPEx ones.
    - Make sure the paths to selection_file, topology, west_home, reference and reference_fop exist and are valid.
    - the reference is the pdb file that will used in EVERY SINGLE progress coordinate calculation.
    - selection_file and reference_fop will be calculated later but still need the names that will be used.
6. Run the westpa_scripts/get_reference_fop.py script. It used the west.cfg as config file. Here is where the selection_file and reference_fop files are generated. 
7. Visually inspect the pocket field of points and the selection string (it uses MDAnalysis syntax) you can re-run the westpa_scripts/get_reference_fop.py script to recalculate them.
8. Change the adaptive_binning/adaptive.py file to reflect the bins minimum and maximum values the progress coordinate must take and the number of bins and dimensions. 
9. Revise env.sh, you need to export the MD engine, this file can get complicated if using a supercomputing center. There are some extra files in the env-crc directory that will help with setup at the CRC (University of Pittsburgh's supercomputing center).
10. Change westpa_scripts/get_pcoord.sh.
    - Make sure to link the files needed to start the simulations.
    - make sure that you are copying all the auxdata files to their respective WESTPA variable.
11. Activate the westpa conda environment and run the init.sh file.
12. If errors occur, check the job_logs directory.
13. Modify the runseg.sh
    - link files necessary files to restart simulation.
    - link all parent auxdata files as parent_AUXDATA.txt
    - make sure the random seed number gets added to the configuration file (comment out the one you are not using).
    - make sure the way to call the MD engine is correct.
    - make sure to call the right trajectory file for the jdistance.py script.
    - make sure that you are copying all the auxdata files to the westpa variables.
14. Modify the MD configuration file in reference (check names and other parameters)
    - Make sure the number of frames corresponds to pcoordlength - 1 (pcoordlength is in the adaptive_binning/adaptive.py file)
15. run the run.sh file (unless you are running it in a supercomputing center).

Notes:
There are a LOT of parts that need to be perfect for it to run, since WESTPA sims are not that easy to set up. To debug stuff check the west.log file and also the output in the job_logs directory.

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
