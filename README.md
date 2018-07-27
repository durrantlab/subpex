# SubPEx code

## What is it?

SubPEx is a tool that uses WESTPA to obtain an ensamble of conformations of the pocket of a protein. This then can be used to do docking .

## How to use it?

This code ONLY works for users with NAMD2.
First you need to load the westpa environment. (e.g. source activate westpa-2017.10)
After thet you will need to run the init.sh script (sometimes you need to run it twice, this needs to be checked)
Once you don't have errors while running the script run the run.sh script. 

For more information on these scripts go to the comments in the files.

## Notes

This a first WORKING test for SubPEx code.

Things that need to be changed:

- Check glitchyness of init.sh 
- Check why MD simulations are "jumping"
- Check why Subpex values are really high
- add more data to the west.h5 file (like pocket volume)
- add a file with all the different pocket shapes
- Add the code to do more than one calculation of pcoordinate per worker

## Important scripts and files and what they do 

__env.sh__ sets up some enviroenmental variables.

__init.sh__ initializes the run. Creates the basis state, the initial states and calculates progress coordiantes for those.

__gen_istate.sh__ copies necessary files to generate initial states.

__get_pcoord.sh__ calls on personall scripts to calculate progress coordinate. This calls for the SubPEx code.

__run.sh__ starts the run.

__runseg.sh__ WESTPA runs this script for each trajectory segment. has three jobs:
- Will create the files necessary for md simulations
- Run the MD simualtion
- Calculate the progress coordiante.

__west.cfg__ the file that contains the configuration of the run.
