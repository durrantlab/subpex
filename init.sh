#!/bin/bash
#
# init.sh
#
# Initialize the WESTPA simulation, creating initial states (istates) and the
# main WESTPA data file, west.h5. 
#
# If you run this script after starting the simulation, the data you generated
# will be erased!
#

source env.sh

# Make sure that WESTPA is not already running.  Running two instances of 
# WESTPA on a single node/machine can cause problems.
# The following line will kill any WESTPA process that is currently running.
pkill -9 -f w_run

# Make sure that seg_logs (log files for each westpa segment), traj_segs (data
# from each trajectory segment), and istates (initial states for starting new
# trajectories) directories exist and are empty. 
rm -rf traj_segs seg_logs istates west.h5 
mkdir -p seg_logs traj_segs istates job_logs

# Copy over the equilibrated conformation from ./prep to bstates/unbound,
# including coordinates, velocities, and box information
rm -rf bstates/first_conformer
mkdir -p bstates/first_conformer

# Define the arguments for the basis states (used for generating initial
# states; in this case we only have one), and target states (used for
# knowing when to recycle trajectories). In this example, we recycle
# trajectories as they reach the bound state; we focus on sampling  
# the binding process (and not the unbinding process).

BSTATE_ARGS="--bstate-file bstates/bstates.txt"
# In this case we are not recycling walkers so no TSTATE is necessary
#TSTATE_ARGS="--tstate bound,1.0"       

# Initialize the simulation, creating the main WESTPA data file (west.h5)
# The "$@" lets us take any arguments that were passed to init.sh at the
# command line and pass them along to w_init.
$WEST_ROOT/bin/w_init \
  $BSTATE_ARGS \
  $TSTATE_ARGS \
  --segs-per-state 4 \
  --work-manager=$WORKMANAGER "$@"

rm -f binbounds.txt