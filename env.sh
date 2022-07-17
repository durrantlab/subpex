#!/bin/bash
#
# env.sh
#
# This script defines environment variables that are used by other shell
# scripts, both when setting up the simulation and when running the simulation.
#

###############################################################################

# STEP 1: If running SubPEx on a system that uses modules, load the modules
# required to run AMBER or NAMD. Here are some examples, though these are not
# likely to work on your specific system.

# module purge
# module load namd/2.12b1-multicore-CUDA

# module load gcc/10.2.0 openmpi/4.1.1
# module load amber/22

###############################################################################

# STEP 2: Set the environment variables for running the AMBER or NAMD
# executable, per your system setup. Here are some examples, with the NAMD
# example uncommented:

export NPROC=12
export NAMDBIN=$(which namd2)
export NAMD="$NAMDBIN +p$NPROC +idlepoll +setcpuaffinity"

# export NPROC=$(nproc)
# SANDER=pmemd.MPI
# export AMBER="mpirun -n $NPROC $SANDER"

###############################################################################

# STEP 3: Check to make sure that the environment variable WEST_ROOT is set.
# Here, the code '[ -z "$WEST_ROOT"]' will return TRUE if WEST_ROOT is not set,
# causing us to enter the if-statement, print an error message, and exit.

if [ -z "$WEST_ROOT" ]; then
  echo "The environment variable WEST_ROOT is not set."
  echo "Try running 'source westpa.sh' from the WESTPA installation directory"
  echo "This is going to cause problems unless you fix it!!!"
fi

###############################################################################

# STEP 4: Explicitly name our simulation root directory. Similar to the
# statement above, we check if the variable is not set. If the variable is not
# set, the we set it.

if [ -z "$WEST_SIM_ROOT" ]; then
  export WEST_SIM_ROOT="$PWD"
fi

export WEST_SIM_ROOT="$PWD"

###############################################################################

# STEP 5: Set the simulation name. Whereas "WEST_SIM_ROOT" gives us the entire
# absolute path to the simulation directory, running the "basename" command will
# give us only the last part of that path (the directory name).

export SIM_NAME=$(basename $WEST_SIM_ROOT)

###############################################################################

# STEP 6: Set the work manager to processes. See
# https://github.com/westpa/westpa/wiki/Running-WESTPA-in-a-multi-node-environment

#export WORKMANAGER="processes"
export WORKMANAGER="serial"
