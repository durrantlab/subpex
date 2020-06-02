# SubPEx code

## What is it?

Subpocket explorer (SubPEx) is a tool that uses weighted ensemble, as implemented in 
WESTPA, to maximizes pocket conformational search with the goal of obtaining an 
ensemble of protein conformations for the use in ensemble docking.

## How to use it?

This code at the moment only works for users with NAMD2. Make sure you have WESTPA 
installed in your computer and that it calls the right python (python2.7). For the 
calculation of the progress coordinate SubPEx needs python3 with MDAnalysis package.
At the moment manual input is needed for the user to be able to adapt it to the protein
of interest. This will be changed soon.
First run env.sh to load all of the WESTPA variables. Then, run init.sh which will 
calculate the basis state progress coordinates and then create the initial states and 
calculate their progress coordinate. Once this 
After that you will need to run the init.sh script (sometimes you need to run it twice, this needs to be checked)
Once you don't have errors while running the script run the run.sh script. 

For more information on these scripts go to the comments in the files.

## Notes

This a WORKING alpha version of SubPEx.

Things that need to be changed:

- add a file with all the different pocket shapes or a directory with protein 
conformations.
- add script to help with setting things up.  

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