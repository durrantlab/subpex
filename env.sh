#!/bin/bash
#
# env.sh
#
# This script defines environment variables that are used by other shell
# scripts, both when setting up the simulation and when running the simulation.
# It will almost certainly need to be modified for your environment, but we
# provide examples below to get you started.

export ENGINE="NAMD"  # NAMD or AMBER
export MODE="MPI"  # MPI, GPU, or MULTITHREAD

module purge

if [ "$ENGINE" == "NAMD" ] ; then

  #############################################################################

  # STEP 1: If running SubPEx on a system that uses modules, load the modules
  # required to run AMBER or NAMD. Here are some examples, though these are not
  # likely to work on your specific system.

  module load namd/2.12b1-multicore-CUDA

  #############################################################################

  # STEP 2: Set the environment variables for running the AMBER or NAMD
  # executable, per your system setup. Here are some examples, with the NAMD
  # example uncommented:

  export NPROC=12
  export NAMDBIN=$(which namd2)
  export NAMD="$NAMDBIN +p$NPROC +idlepoll +setcpuaffinity"

elif [ "$ENGINE" == "AMBER" ] ; then

  ###############################################################################

  # STEP 1: If running SubPEx on a system that uses modules, load the modules
  # required to run AMBER or NAMD. Here are some examples, though these are not
  # likely to work on your specific system.

  # AMBER
  module load gcc/8.2.0 openmpi/4.0.3
  module load amber/22

  #############################################################################

  # STEP 2: Set the environment variables for running the AMBER or NAMD
  # executable, per your system setup. Here are some examples, with the NAMD
  # example uncommented:

  if [ "$MODE" == "MPI" ] ; then
    # AMBER, MPI
    export NPROC=$(nproc)
    SANDER=pmemd.MPI
    export AMBER="mpirun -n $NPROC $SANDER"
  elif [ "$MODE" == "GPU" ] ; then
    # AMBER, GPU
    SANDER=pmemd.cuda
    export AMBER="$SANDER"
  fi
fi

###############################################################################

# STEP 3: Check to make sure that the environment variable WEST_ROOT is set.
# Here, the code '[ -z "$WEST_ROOT"]' will return TRUE if WEST_ROOT is not set,
# causing us to enter the if-statement, print an error message, and exit.

# Note: Certain supercomputing centers may require you to specify the full path
# to your westpa environment.
echo "Using conda: $(which conda)"
source $(conda info | grep -i 'base environment' | awk '{print $4}')/etc/profile.d/conda.sh
conda activate westpa

if [ -z "$WEST_ROOT" ]; then
  echo "The environment variable WEST_ROOT is not set."
  # echo "Try running 'source westpa.sh' from the WESTPA installation directory"
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
# https://westpa.github.io/westpa/users_guide/wwmgr.html

if [ "$MODE" == "MPI" ] ; then
  export WORKMANAGER="zmq"  # Good for slurm, multiple nodes
elif [ "$MODE" == "GPU" ] ; then
  export WORKMANAGER="serial"  # Good for running on single machine (single GPU).
elif [ "$MODE" == "MULTITHREAD" ] ; then
  export WORKMANAGER="processes"  # Single node with multiple processors. TODO: Not tested.
fi
